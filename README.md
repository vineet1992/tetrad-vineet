### Vineet's PhD Thesis All Code
This is a local version of the tetrad project that builds upon graphical structure learning in the case of mixed data.

For the full tetrad project, please visit: https://github.com/cmu-phil/tetrad

Included in here are a suite of java files for running causal discovery methods that handle latent variables and mixed continuous and discrete datasets. My thesis is focused upon incorporating prior knowledge into these methods for biomedical research applications. Many of the outputs of this project are focused on analzying observational, static, multi-omics datasets. 

Here is a pictoral description of my thesis:

![Image of Thesis](https://github.com/vineet1992/tetrad-vineet/blob/master/Thesis%20Workflow.png)

# Steps to use the releases of this thesis

### 1. Install Java 

In order to run these methods you will need to have Java installed on your machine, please see these instructions for doing so: https://www.informationweek.com/desktop/how-to-install-java-runtime-environment-in-windows/d/d-id/1099686

### 2. Set Java Home

Windows 10 and Windows 8
In Search, search for and then select: System (Control Panel)
Click the Advanced system settings link.
Click Environment Variables. In the section System Variables, find the PATH environment variable and select it. Click Edit. If the PATH environment variable does not exist, click New.
In the Edit System Variable (or New System Variable) window, specify the value of the PATH environment variable. Click OK. Close all remaining windows by clicking OK.
Add the following to the PATH variable (**make sure this is separated by a semicolon from everything that was already in there**)
"C:/Program Files/Java/jre1.8.0_161/bin;" 

If you have a different version of Java then you need to change the 1.8.0_161 to your version, navigate to this directory on your computer to double check!

For Mac OS, you usually don't need to do anything


Next be absolutely sure that your .jar file, your datasets, and all files that will be used by your .jar file are in the same folder on you computer. 

Then, navigate to this directory on your windows command prompt or mac terminal by using the command "cd" 
Example: cd C:/Users/vineet/Documents/Directory_With_My_Files





