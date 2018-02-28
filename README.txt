Description of solution by:
Marek Cygan
handle: marek.cygan
email: marek.cygan@gmail.com

Please follow the instructions when running the scripts:
a) mount both ssd drives to some directories, say ~/ssd1 and ~/ssd2,
b) copy all the files from the /home/ubuntu/mydata to _BOTH_ ssd directories,
c) go to the connectivity-map2 directory:
cd /home/ubuntu/connectivity-map2
d) make clean
e) make
c) use the script:
./run.sh upfile downfile outputfile ssd1_directory ssd2_directory
Important note: in the command above please specify the output file to be on the same ssd drive as ssd1_directory is.

Input/output assumptions:
a) the output file is written as a list of fp32 numbers in the same order as in the problem statement,
b) the up and down files are assumed to contain the name column as in the ones provided in the problem statement.
