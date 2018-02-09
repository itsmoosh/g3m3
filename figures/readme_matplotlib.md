Below includes the steps to produce an eps figure from cut multifluid data. 

    0 - Ensure you have the following installed:
        
        python3
        python3-matplotlib
        python3-numpy
        python3-scipy

    1 - Ensure your data is cut to an ASCII file in the format:

            (( example files included in ./data/ directory ))
            
            X1   Y1   var1    var2    var3    ...   varN
            X2   Y2   var1    var2    var3    ...   varN
            ... etc.

    2 - From the command line, execute:

        python3 alfven_example.py

    3 - To open the file, execute:

        evince plot_alfven_example.eps &

            (( Note: with '&' specified you can keep the image open and iterate script changes ))
