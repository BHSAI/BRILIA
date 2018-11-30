# BRILIA file and web access details

When using BRILIA, you make be prompted to accept permission for BRILIA.exe to access/modify files. Pushing "cancel" can work, but to better inform the IT security team, here is a list file or website the BRILIA may access.

## Files that are added:  

  * Databases folder containing the VDJ sequences  
  * Examples folder that contains examples sequences for BRILIA to test  
  * Any output files for your data  
  * Temp files when using the EXE, which can be found via:  
     ```  
     BRILIA> findRoot print  
     ```  

## Files that can be removed or modified:  

  * Any output files that are overridden when redoing the analysis

## Website that may be accessed:  

  * BRILIA can check https:/github.com/BHSAI/BRILIA/README.md to see what version is the latest, but only if invoking the following command:  
    ```  
    BRILIA> checkLatestVersion  
    ```  