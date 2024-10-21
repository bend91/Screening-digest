# Screening-digest
Python tool for identifying optimal enzymes to use for a screening digest.

run "python functions.py" from command line, you will be prompted for two DNA sequences, either paths to snapgene files or manually enter the sequences. The program will find enzymes that cut in an optimal way that will distinguish between the two sequences. Currently assumes all sequences are plasmids, will introduce functionality to choose either circular or linear sequences.
With manual entry, if copying and pasting ensure the terminal accepts more characters to be pasted than your sequence length, I ran into issues with testing as my terminal (Terminus add on in sublime text) would only accept 1024 characters being pasted which is obviously a lot smaller than a plasmid sequence.

##GUI 
python main.py 

A work in progress using Pyside, enter the two sequences in the text box, will give a list of enzymes and show the digest simulation when you click on one.

Future developments
- Make an executable for the GUI so python and pyside aren't required
- make the GUI accept drag and drop of snapgene files
- Allow saving of the gel simulation in the GUI
- Implement more than one comparison and include enzymes that cut in one sequence and not the other
- Allow a custom list of enzymes to check
