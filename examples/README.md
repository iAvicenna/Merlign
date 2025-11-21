# Mer|ign Script

This folder contains a bin script called merlign which is a fully fledged 
bash script that can accept its own inputs and has its own internal tests. 
For more information on parameters call

    ./merlign -h

For more information on how to use the script and how to structure the project
folder call

    ./merlign -i
    
This script can use a settings file, an example of which is given in merlign_settings.txt.
An example usage would be:

    #!/bin/bash

    name="myproject"
    seqs_folder="/home/seqs_folder/"
    proj_dir="./"

    merlign -s merlign_settings.txt -Q 30 -N 6 "${seqs_folder}" "${proj_dir}" "$name" > "$name.merlign.log"

There is an accessory library called [mitril](https://github.com/iAvicenna/mitril) for analysing the results 
produced by Merlign. 
