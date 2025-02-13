# XML data project

> AIM 1: Write function to generate fasta including sequences of all variants from UniProt variant XML files.

----------------------------

### November 5, 2024

I looked into the structure of XML files, and researched about possibilities to read and write XML files using Python.

The XML files for protein entries from UniProt does not seem to include the variant data included in the variant XML, as the categories such as "variantlocation" are not found in the protein entry XML.

----------------------------

### November 21, 2024

I looked into ways to handle XML files using Python and found ElementX. I looked into how it works. I looked into how to save everything into fasta files the easiest. I am going to use BioPython for that.

I also looked into the different types of variants there are. I found missense, ...

I wrote a pseudocode for AIM1.

----------------------------

### November 22, 2024

Implemented the first few lines of pseudocode.

----------------------------

### November 25, 2024

Continued implementing pseudocode. Stepped on issue with for loop not correctly iterating through features apparently (features do not have children they should have). 

----------------------------

### November 26, 2024

Fixed issue of iteration not working properly. The mistake was using find instead of findall to access features.

----------------------------

### November 29, 2024

I finished writing a first version of the function to create a variant fasta from the variant XML (from UniProt). The function is called `fasta_from_variantsXML()` (in script `XML.py`). The function is designed to work on all UniProt variants XML files by simply providing the path to the XML file (and optionally the CDS nucleotide sequence of the protein), and writes a fasta file including all variant sequences as output. I tested and debugged the function (only) on the (UniProt) XML of interest including all variants of the human transthyretin (accession: P02766). The function could successfully create a fasta file with (almost) all variants included in the XML (also the right sequences after debugging, e.g. recognizing falsely as missense annotated consequence of variant, e.g. handling variants differing in multiple positions correctly). Variants with consequence types "frameshift" (thus CDS nucleotide sequence can be provided but processing frameshifts not implemented yet) and "stop lost" (not inferable what comes after stop codon, perhaps partially possible in future by allowing to provide whole gene sequence). Also, the script most likely cannot recognize all possible consequence types yet and further types should be added when encountered in future.

I will send the transthyretin variants fasta to Christian. All variants except one framshift variant (position 91) and one stop lost variant are included (see reasons above).

----------------------------

### December 2, 2024

I talked to Christian about the file I sent and the general project. The script could be improved in the following ways:
- the original sequence of the protein should be included at position 1 in the generated fasta
- more consequenceTypes should be added based on tests done with other XML files
- improve user-friendliness (so that the script is easily usable by all people in the group)
- advanced: implement functionality to iterate through database of many different proteins, download associated variants XML automatically, and create variant fasta for each

Today, I already added that now the original sequence will also be in the generated fasta (with right database at start, sp or tr, based on review argument). I also added the possibility to provide an output directory as argument.

Further, I tested some other XMLs to search for more consequenceTypes or weird behavior. Based on one file, I added the abiliy to recognize defect start codons (variant is skipped). I now downloaded the longest protein (curated) there is in humans (Q8WZ42, Titin) and ran the function for it. Many "Unknown consequenceType" messages were printed. I can go through each one of them (by breakpoint) on the next day to add the ability to recognize these missing consequenceTypes.

----------------------------

### December 5, 2024

After testing several different variant XMLs from UniProt entries, I implemented the recognition and handling of two more consequenceTypes:
- **"-"** (apparently translates to delins variants): indicated amino acid(s) are deleted and the indicated new amino acids are inserted at the same postion
- **"stop retained"**: no change in amino acid sequence occurs (only prints message that variant is present)

For now, I prevented the message "Frameshift variant cannot be processed. No CDS provided." from being printed to avoid confusion of the function being able to handle frameshifts in the first place. Instead, a new message appears that just informs about a framshift variant at the respective position.

----------------------------

### December 17, 2024

***Aim: Make script more user-friendly.***

I started implementing a GUI using the Python-included package **Tkinter**. This package allows for a simple, easy to use GUI to select files and submit, and is also functioning across platforms. I set up the template of the GUI implementation and will continue modifying it to the needs of the script on the next day.

----------------------------

### December 18, 2024

I modified and added to the Tkinter implementation based on the needs of the script. For example, I modified the displayed names, wanted input file type, and set a standard path for the output directory being the current working directory. This concludes the work on user-friendliness for now (if needed, more will be added).

Next, I should work on an implementation to accept normal XML files (not variant XML files, and also allowing multiple entries), and automatic search for each variant XML file to process and generate a fasta.

----------------------------

### February 12, 2025

***Aim: Implement accepting normal XML files and running for all entries***

- added recognition of XML file type
- added download of variant XML file function
- implemented allowing multiple file selection (in Tkinter) and processing (works with variant XML files)
- started implementing processing of normal XML files in general

Next time:
- finish normal XML files handling
- fix Tkinter window not closing when pressing submit
- fix only allowing one time browse of input files in Tkinter

----------------------------

### February 13, 2025

***Aim: Release first version able to handle multiple and normal UniProt XML files***

- finished implementation of multiple and normal XML file handling
- modified Tkinter GUI for user-friendliness:
    - fixed issue of one time browse of input files
    - added console window and redirected print commands to it in real-time using threading
    - added close button
- modified and added several print messages

This concludes the first version of FASTA2XML.