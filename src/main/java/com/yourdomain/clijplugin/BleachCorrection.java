package com.yourdomain.clijplugin;

/* TODO:
    -fix exponential method
    	-returns images that are all maximum intensity
    -add histogram method
    -convert to cli2
    -handle multichannel and timeseries data
*/

/*
 * Bleach correction on the GPU
 *
 * @Author: Robert Haase, rhaase@mpi-cbg.de
 * July 2019
 *
 * Code partly adapted from:
 *
 * Bleach Correction with three different methods, for 2D and 3D time series.
 * 	contact: Kota Miura (miura@embl.de)
 *
 * 	Simple Ratio Method:
 * 		Plugin version of Jens Rietdorf's macro, additionally with 3D time series
 * 			see http://www.embl.de/eamnet/html/bleach_correction.html
 *
 *  Exponential Fitting Method:
 *  	Similar to MBF-ImageJ method, additionally with 3D time series.
 *  		See http://www.macbiophotonics.ca/imagej/t.htm#t_bleach
 *  	MBF-ImageJ suggests to use "Exponential" equation for fitting,
 *  	whereas this plugin uses "Exponential with Offset"
 *
 *  HIstogram Matching Method:
 *  	This method does much better restoration of bleaching sequence
 *  	for segmentation but might not good for intensity quantification.
 *  	See documentation at http://cmci.embl.de
 *
 *
 * Copyright © 2010 Kota Miura
 * License: GPL 2
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 2
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import emblcmci.BleachCorrection_ExpoFit;
import emblcmci.BleachCorrection_MH;
import emblcmci.BleachCorrection_SimpleRatio;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.measure.CurveFitter;
import ij.plugin.Duplicator;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.Fitter;
import ij.process.ImageProcessor;

// added CLIJ imports
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij.macro.AbstractCLIJPlugin;
import net.haesleinhuepf.clij.macro.CLIJMacroPlugin;
import net.haesleinhuepf.clij.macro.CLIJOpenCLProcessor;
import net.haesleinhuepf.clij.macro.documentation.OffersDocumentation;
import org.scijava.plugin.Plugin;

import java.util.Arrays;
import java.util.HashMap;

@Plugin(type = CLIJMacroPlugin.class, name = "CLIJ_bleachCorrection")
public class BleachCorrection extends AbstractCLIJPlugin implements CLIJMacroPlugin, CLIJOpenCLProcessor, OffersDocumentation {

	String[] correctionMethods =  { "Simple Ratio", "Exponential Fit"}; //, "Histogram Matching" };

	@Override
	public boolean executeCL() {
		// get input from user; a dialog is opened before this method is called
		ClearCLBuffer input = (ClearCLBuffer)( args[0]);
		ClearCLBuffer output = (ClearCLBuffer)( args[1]);
		String method = args[2].toString().toLowerCase();

		// go through available methods and search for the entered one
		int selectedMethod = 0;
		for (int i = 0; i < correctionMethods.length; i++) {
			if (correctionMethods[i].toLowerCase().compareTo(method) == 0) {
				selectedMethod = i;
				break;
			}
		}

		// measure SUM intensities for each slice
		double[] sumIntensities = clij.op().sumPixelsSliceBySlice(input);
		float[] correctionFactors = new float[sumIntensities.length];

		// determine correction factors
		if (selectedMethod == 0) { // simple ratio
			// adapted from https://github.com/fiji/CorrectBleach/blob/master/src/main/java/emblcmci/BleachCorrection_SimpleRatio.java#L109
			double referenceInt = sumIntensities[0];
			for (int i = 0; i < input.getDepth(); i++) {
				double currentInt = sumIntensities[i];
				float ratio = (float)(referenceInt / currentInt);
				correctionFactors[i] = ratio;
				IJ.log("frame"+ Integer.toString(i+1) + "mean int="+ currentInt +  " ratio=" + ratio);
			}
		} else if (selectedMethod == 1){ // exponential fit
			double[] yA = sumIntensities;
			double[] xA = new double[yA.length];
			long pixelsPerPlane = (input.getWidth() * input.getHeight());
			System.out.println("using int");

			for (int i = 1; i < input.getDepth(); i++) {
				yA[i] = yA[i] / pixelsPerPlane;
				xA[i] = i;
			}

			// adapted from https://github.com/fiji/CorrectBleach/blob/master/src/main/java/emblcmci/BleachCorrection_ExpoFit.java#L72
			CurveFitter cf = new CurveFitter(xA, yA);
			double firstframeint = yA[0];
			double lastframeint = yA[yA.length-1];
			double guess_a = firstframeint - lastframeint;
			if (guess_a <= 0){
				IJ.error("This sequence seems to be not decaying");
				return false;
			}
			double guess_c = lastframeint;
			double maxiteration = 2000;
			double NumRestarts = 2;
			double errotTol = 10;
			double[] fitparam = {-1*guess_a, -0.0001, guess_c, maxiteration, NumRestarts, errotTol};

			cf.setInitialParameters(fitparam);
			cf.doFit(11); //
			Fitter.plot(cf);
			IJ.log(cf.getResultString());

			double[] respara = cf.getParams();
			System.out.println(respara);
			double res_a = respara[0];
			double res_b = respara[1];
			double res_c = respara[2];

			for (int i = 0; i < input.getDepth(); i++){
				correctionFactors[i] = (float)(calcExponentialOffset(res_a, res_b, res_c, 0.0) / calcExponentialOffset(res_a, res_b, res_c, (double) (i + 1)));
			}
		} else {
			IJ.log("Method " + selectedMethod + "not supported yet.");
			return true;
		}

		return clij.op().multiplySliceBySliceWithScalars(input, output, correctionFactors);
	}

	/** calculate estimated value from fitted "Exponential with Offset" equation
	 *
	 * @param a  magnitude (difference between max and min of curve)
	 * @param b  exponent, defines degree of decay
	 * @param c  offset.
	 * @param x  timepoints (or time frame number)
	 * @return estimate of intensity at x
	 */
	public double calcExponentialOffset(double a, double b, double c, double x){
		return (a * Math.exp(-b*x) + c);
	}

	@Override
	public String getParameterHelpText() {
		return "Image source, Image destination, String method";
	}

	@Override
	public String getDescription() {
		return "Bleach correction; translated from Kota Muiras version: \n" +
				"https://github.com/fiji/CorrectBleach\n\n" +
				"Enter one of these methods: " + Arrays.toString(correctionMethods);
	}

	@Override
	public String getAvailableForDimensions() {
		return "2D, 3D";
	}
}
