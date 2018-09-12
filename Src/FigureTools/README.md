# FigureTools v1.1
MATLAB codes for common plotting tasks, such as saving at 600 DPI or vector graphics, inverting image color for intramural seminar, adding panel labels to subplots, etc.

### Matlab Usage
  1. Download all contents
  2. In Matlab, add the path to FigureTools folder via `addpath('path_to_folder\FigureTools')`
  3. To see if it works, run the following commands
  ```
    Gx = loadDemoPlot;
    savePlot(Gx, 'SaveAs', '')
  ```
  NOTE: MS Office no longer supports .eps vector graphics due to security vulnerability issues. You must use .emf files instead.

### What do you want do?

  * **See what each code does?**  
    [`showDemo`](https://github.bhsai.net/dlee/FigureTools/blob/master/showDemo.m)
    
  * **Save figures according a set DPI and file format?**   
    [`savePlot`](https://github.bhsai.net/dlee/FigureTools/blob/master/savePlot.m)

  * **Invert color of image files? Good for INTRAMURAL SEMINARS!**  
    [`invertPicColor`](https://github.bhsai.net/dlee/FigureTools/blob/master/invertPicColor.m)
    
  * **Invert color of a matlab figure?**  
    [`invertFigColor`](https://github.bhsai.net/dlee/FigureTools/blob/master/invertFigColor.m)

  * **Remove excess white space borders and/or rescale subplots in a figure?**  
    [`resizeSubplots`](https://github.bhsai.net/dlee/FigureTools/blob/master/resizeSubplots.m)
    
  * **Set the number of decimals in the X or Y axes?**  
    [`setPlotTickDecimal`](https://github.bhsai.net/dlee/FigureTools/blob/master/setPlotTickDecimal.m)
    
  * **Add labels "a)", "b)", ... on subplots?**  
    [`labelSublots`](https://github.bhsai.net/dlee/FigureTools/blob/master/labelSubplots.m)
    
  * **Resize the figure to journal specificiations?**  
    [`resizeFigure`](https://github.bhsai.net/dlee/FigureTools/blob/master/resizeFigure.m)
    
  * **Reposition figures going outside your monitor top to the center of your monitor?**   
    [`centerFigureOnMonitor`](https://github.bhsai.net/dlee/FigureTools/blob/master/centerFigureOnMonitor.m)

  * **Get the current figure handle without creating a new figure by accident?**  
    [`getOpenGCF`](https://github.bhsai.net/dlee/FigureTools/blob/master/getOpenGCF.m)
    
  * **Get the current axes handle without creating a new axes by accident?**  
    [`getOpenGCA`](https://github.bhsai.net/dlee/FigureTools/blob/master/getOpenGCA.m)

  * **Set a property of ALL AXES in a figure?**  
    [`setAxes`](https://github.bhsai.net/dlee/FigureTools/blob/master/setAxes.m)
