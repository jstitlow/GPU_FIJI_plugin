input = getTitle();
getDimensions(width, height, channels, slices, frames);
output = input+"zcorr";; 

// Init GPU
run("CLIJ Macro Extensions", "cl_device=");
Ext.CLIJ_clear();

// push images to GPU
Ext.CLIJ_push(input);

// cleanup ImageJ
//run("Close All");

// Blur in GPU
//Ext.CLIJ2_gaussianBlur3D(input, blurred, 10, 10, 1);

// Z correction in GPU
Ext.CLIJ_bleachCorrection(input, output, "Simple Ratio");

// Get results back from GPU
Ext.CLIJ_pull(output);

// Cleanup by the end
Ext.CLIJ_clear();