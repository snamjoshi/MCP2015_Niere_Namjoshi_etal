---
title: "ImageJ Macros"
author: "Sanjeev V Namjoshi"
date: "Revised for GitHub: October 16, 2015"
output: 
  html_document: 
    highlight: haddock
    keep_md: yes
---

**Macros written in late July and early August 2015**  

The following is a collection of ImageJ macros I used in some of the analysis sections of the paper. Many of these macros will need tweaking to make them work for your purposes. 

SEE ALSO: JacopAnalysis.rmd for code to tidy up the output from the ImageJ log after running JACoP.

### Calculate RGB values for all images in a specified directory

    dir = getDirectory("Choose a Directory");
    list = getFileList(dir);
    
    for(i = 0; i < list.length; i++) {
	    open(dir+list[i]);
          print(getTitle());
          run("RGB Measure");
	    run("Close All");
    }

### Threshold RGB values to specified amount

The first parameter in the "setMinAndMax" function below specifies the intensity of the channel.

    dir=getDirectory("Choose a Directory");
    print(dir);
    thresholded=dir + "\\thresholded\\";
    print(thresholded);
    File.makeDirectory(thresholded);
    list = getFileList(dir);

    for (i=0; i<list.length; i++) {
       open(dir+list[i]);
       imgName = getTitle();
       baseNameEnd=indexOf(imgName, ".tif");
       baseName=substring(imgName, 0, baseNameEnd);
       setMinAndMax(0, 255, 4);  //Red Channel
       setMinAndMax(41, 255, 2);  //Green Channel
       setMinAndMax(20, 255, 1);  //Blue Channel
       selectWindow(imgName);
       save(thresholded+imgName);
       run("Close All"); 
    }
	
### Automate JACoP analysis

"Thra" and "Thrb" indicates the threshold for image A and image B respectively. The thresholding first acts on the red/blue channel, the red/green channel, and finally the green/blue channel. These orders can be altered to suit your needs.

    dir=getDirectory("Choose a Directory");
    list = getFileList(dir);
    Array.sort(list);
     
    for (i=0; i<list.length; i+=3) {
	    open(list[i]);
	    original=getTitle;
      blue=substring(original,0,lengthOf(original)-9)+"-blue.tif";     
    
      open(list[i+1]);  
	    original=getTitle;
      green=substring(original,0,lengthOf(original)-10)+"-green.tif";

      open(list[i+2]);
      original=getTitle;
      red=substring(original,0,lengthOf(original)-8)+"-red.tif";
    
      run("JACoP ", "imga="+red+" imgb="+blue+" thra=648 thrb=517 pearson mm");
      run("JACoP ", "imga="+red+" imgb="+green+" thra=648 thrb=517 pearson mm");
      run("JACoP ", "imga="+green+" imgb="+blue+" thra=648 thrb=517 pearson mm");

      close(red);
      close(green);
      close(blue);
    }

### Quick Rotate

A simple macro to automatically rotate the image at a specified angle. Saves time from having to open the menu and select the options.

    run("Rotate... ", "angle=-20 grid=1 interpolation=None");

### Split images into multiple channels and save each image

    dir=getDirectory("Choose a Directory");
    print(dir);
    splitDir=dir + "\Split\\";
    print(splitDir);
    File.makeDirectory(splitDir);
    list = getFileList(dir);

    for (i=0; i<list.length; i++) {
         if (endsWith(list[i], ".tif")){
           print(i + ": " + dir+list[i]);
           open(dir+list[i]);
           imgName=getTitle();
           baseNameEnd=indexOf(imgName, ".tif");
           baseName=substring(imgName, 0, baseNameEnd);

           run("Split Channels");
           selectWindow("C3-" + imgName);
           rename(baseName + "-blue.tiff");
           saveAs("Tiff", splitDir+baseName + "-blue.tif");
           close();
           selectWindow("C2-" + imgName);
           saveAs("Tiff", splitDir+baseName + "-green.tif");
           close();
           selectWindow("C1-" + imgName);
           saveAs("Tiff", splitDir+baseName + "-red.tif");
	         close();
	         selectWindow("C4-" + imgName);
           saveAs("Tiff", splitDir+baseName + "-gray.tif");

           run("Close All");
        }
    } 

### Multi-Channel Measure

Measures RGB values for all channels in a directory (not just RGB)

    dir = getDirectory("Choose a Directory");
    list = getFileList(dir);

    for(i = 0; i < list.length; i++) {
	    open(dir+list[i]);
      print(getTitle());
      run("Measure");
	    run("Close All");
    }










