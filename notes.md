# BleachCorrection plugin- notes  
## mvn install
* required JDK1.8
	* made an alias to call mvn with JDK1.8
	* alias mvn1='JAVA_HOME=`/usr/libexec/java_home -v 1.8` && mvn'
* issue with duplicate classes 
	* found the duplicate classes in External Libraries in the IntelliJ project
	* workaround by mvn install -Denforcer.skip=true
* test macro in FIJI returns "Fouhnd 1 method..."!!!!
