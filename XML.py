import xml.etree.ElementTree as ET
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np
import sys
import requests
import threading

# ------------------------------------------


def download_variantXML(accession, outdir, max_acc, i):

    url = f'https://www.ebi.ac.uk/proteins/api/variation/{accession}?format=xml'
    outpath = os.path.join(outdir, f'{accession}.xml')

    try:
        # make download request
        resp = requests.get(url, timeout= 10)
        resp.raise_for_status()

        # write xml
        with open(outpath, 'wb') as file:
            file.write(resp.content)
        if i == max_acc:
            print(f'({i}/{max_acc}) Successfully downloaded {accession}.xml\n', flush= True)
        else:
            print(f'({i}/{max_acc}) Successfully downloaded {accession}.xml', flush= True)

    except requests.exceptions.RequestException as e: 
        if i == max_acc:
            print(f'{i}/{max_acc}) Failed to download {accession}.xml: {e}\n', flush= True)
        else:
            print(f'{i}/{max_acc}) Failed to download {accession}.xml: {e}', flush= True)

### Pseudocode (main analysis per entry)
# Read XML file of TRR and set it up
# Save original sequence
# Find all features
# iterate through all features
    # find position and save
    # find variant
    # find consequencetype and save
    # decide what to do based on consequence type
        # missense:
            # save wildtype and mutatedtype
            # replace amino acid in sequence
        # stop gained:
            # cut sequence at current position
        # insertion
            # check variantLocation
                # dup inside: insert what is in wildtype another time
                # no dup: print("Unknown insertion. Skipped")
        # in-frame deletion
            # remove current amino acid
        # frameshift
            # if nucleotide sequence provided:
                # calculate new sequence from this position until stop
                # if nop stop codon:
                    # print("Frameshift leads to no stop. Skipped.")
            # else:
                # print("Frameshift, no nucleotide sequence provided")
        # stop lost
            # print("Stop list variant. Skipped.")
    # save everything into list with metadata
    # save as fasta


def fasta_from_variantsXML(files, outdir = '.', CDS = None, reviewed = True):

    # evidence level dict (UniProt)
    ev_dict = {'Evidence at protein level': 1, 'Evidence at transcript level': 2, 'Inferred from homology': 3, 'Predicted': 4, 'Uncertain': 5}

    ### multiple file handling (files saves paths) -> create fpaths list of variant xml files for processing in next loop
    # entire list of paths to process
    fpaths = []
    for filepath in files:

        # read tree
        tree = ET.parse(filepath)
        # find root
        root = tree.getroot()

        # check what is xml format
        if root.tag.endswith('uniprot'):

            print('Uniprot XML detected. Downloading associated variant XML files...', flush= True)

            namespace = {"ns": root.tag.split("}")[0].strip("{")}
            entries = root.findall('.//ns:entry', namespace)
            
            # create accession list for download
            accessions = []
            for entry in entries:
                acc = entry.find('ns:accession', namespace).text
                accessions.append(acc)

            max_acc = len(accessions)

            # make folder for downloads
            download_path = f'{outdir}/downloads'
            if not os.path.exists(download_path):
                os.mkdir(download_path)
            
            # download variant xmls
            for i, accession in enumerate(accessions):
                download_variantXML(accession, download_path, max_acc, i+1)
                fpaths.append(f'{download_path}/{accession}.xml')


        elif root.tag.endswith('entryFeature'):
            # append raw filepath to process (already variant xml file)
            fpaths.append(filepath)
        
        else:
            print(f'Could not process "{filepath}"... Wrong xml file type.', flush= True)
            continue
    

    print('------------------------------------------\nStarting XML2FASTA conversions...\n------------------------------------------', flush= True)

    f_max = len(fpaths)

    # process
    for f, filepath in enumerate(fpaths):
        
        print(f'({f+1}/{f_max}) Processing {os.path.basename(filepath)}...', flush= True)

        # prepare sequence/fasta dict
        fasta_seqs = []

        # read tree
        tree = ET.parse(filepath)
        # find root
        root = tree.getroot()

        # namespace
        namespace = {"ns": root.tag.split("}")[0].strip("{")}
        
        # save metadata
        name = root.find('ns:name', namespace).text
        accession = root.find('ns:accession', namespace).text
        protname = root.find('ns:proteinName', namespace).text
        genename = root.find('ns:geneName', namespace).text
        protein_existence = ev_dict[root.find('ns:proteinExistence', namespace).text]
        sequence = root.find('ns:sequence', namespace).text
        sequence_version = root.find('ns:sequence', namespace).get('version')
        taxid = root.find('ns:taxid', namespace).text
        organism = root.find('ns:organismName', namespace).text

        # database based on review status
        if reviewed == True:
            dbase = 'sp'
        else:
            dbase = 'tr'

        # add original sequence in beginning
        fasta_seqs.append(SeqRecord(Seq(sequence), id= f'{dbase}|{accession}|{name}', description= f'{protname} OS={organism} OX={taxid} GN={genename} PE={protein_existence} SV={sequence_version}'))

        # find all features (one feature per variant), later one sequence per feature
        features = root.findall('.//ns:feature', namespace)

        # iterate through features
        for feature in features:
            
            # find position of current feature
            single_position = 0
            position_feat = feature.find('.//ns:position', namespace)
            if position_feat is not None:
                position = int(position_feat.get('position'))
                single_position = 1
            else:
                position = int(feature.find('.//ns:begin', namespace).get('position'))
                position2 = int(feature.find('.//ns:end', namespace).get('position'))

            # find variant
            variant = feature.find('.//ns:variant', namespace)

            # get consequenceType inside variant
            conseq = variant.find('ns:consequenceType', namespace).text

            # access wildType
            wildtype = variant.find('ns:wildType', namespace).text

            # access mutatedType (sometimes not available, then None)
            mutated = variant.find('ns:mutatedType', namespace)

            # define variable for conditions to avoid false processing of wrongfly annotated variants (as missense)
            stop_gained = 0
            if mutated is not None:
                if mutated.text == '*':
                    stop_gained = 1
                elif mutated.text == '?' and position == 1:
                    print(f'Variant at position {position}: Start codon defective. Skipped.', flush= True)
                    continue

            # decide what to do based on consequenceType
            if conseq == 'missense' and stop_gained == 0: # some inframe deletions are wrongfly annotated as missense (extra variable stop_gained to check)
                # replace amino acid at position and save new sequence
                seq_new = list(sequence)
                seq_new[position-1] = mutated.text
                # define description for fasta header
                desc = f'{wildtype}{position}{mutated.text}'
            elif conseq == 'stop gained' or stop_gained == 1: # extra variabel stop_gained to detect these variants when wrongly annotated
                # define new sequence from start to current position
                seq_new = list(sequence)
                seq_new = seq_new[:position-1]
                # define description for fasta header
                desc = f'{wildtype}{position}*'
            elif conseq == 'insertion':
                # check variantLocation
                varlocs = variant.findall('.//ns:variantLocation', namespace)
                # iterate through varloc
                for varloc in varlocs:
                    loc = varloc.get('loc')
                    if 'dup' in loc:
                        if single_position == 1:
                            # insert wildtype again
                            seq_new = list(sequence)
                            seq_new.insert(position, wildtype)
                            # define description for fasta header
                            desc = f'{wildtype}{position}dup'
                        else:
                            # insert all wildtypes again
                            wildtypes = list(wildtype)
                            seq_new = list(sequence)
                            desc = ''
                            for i, wt in enumerate(wildtypes):
                                seq_new.insert(position+i+len(wildtypes)-1, wt)
                                # define description for fasta header
                                desc = desc + f'{wt}{position+i}_'
                            desc = desc[:-1] + 'dup'
                        break
                    else:
                        continue
                # if no new sequence was generated
                if not seq_new:
                    print(f'Variant at position {position}: Unknown insertion. Skipped.', flush= True)
                    continue
            elif conseq == 'inframe deletion':
                if single_position == 1:
                    # remove current amino acid
                    seq_new = list(sequence)
                    seq_new.pop(position-1)
                    # define description for fasta header
                    desc = f'{wildtype}{position}del'
                else:
                    # remove all amino acids
                    wildtypes = list(wildtype)
                    seq_new = list(sequence)
                    desc = ''
                    for i, pos in enumerate(np.arange(position, position2+1)):
                        seq_new.pop(pos-1-i) # -i because sequence get 1 shorter with every iteration
                        # define description for fasta header
                        desc = desc + f'{wildtypes[i]}{pos}_'
                    desc = desc[:-1] + 'del'
            elif conseq == 'frameshift':
                # check if nucleotide sequence provided (CDS)
                print(f'Variant at position {position}: Frameshift variant. Skipped', flush= True)
                continue
                if CDS:
                    print(f'Variant at position {position}: Frameshifts not implemented yet. Skipped.')
                    continue
                else:
                    print(f'Variant at position {position}: Frameshift variant cannot be processed. No CDS provided. Skipped.')  
                    continue      
            elif conseq == '-':
                # check variantLocation
                varlocs = variant.findall('.//ns:variantLocation', namespace)
                # iterate through varloc
                for varloc in varlocs:
                    loc = varloc.get('loc')
                    if 'delin' in loc:
                        if single_position == 1:
                            # delete wildtype, insert mutatedType starting at same position
                            seq_new = list(sequence)
                            seq_new.pop(position-1)
                            mutated_ins = list(mutated.text)
                            for i, mt in enumerate(mutated_ins):
                                seq_new.insert(position+i-1, mt)
                            # define description for fasta header
                            desc = f'{wildtype}{position}delins{mutated.text}'
                        else:
                            # delete all wildtype, insert mutatedType starting at beginning position
                            seq_new = list(sequence)
                            del seq_new[position-1:position2]
                            mutated_ins = list(mutated.text)
                            for i, mt in enumerate(mutated_ins):
                                seq_new.insert(position+i-1, mt)
                            desc = f'{wildtype[0]}{position}_{wildtype[-1]}{position2}delins{mutated.text}'
                        break
                    else:
                        continue
                # if no new sequence was generated
                if not seq_new:
                    print(f'Variant at position {position}: Type of special variant unknown (consequenceType = "-"). Skipped.', flush= True)
                    continue
            elif conseq == 'stop lost':
                print(f'Variant at position {position}: Stop lost variant. Skipped.', flush= True)
                continue
            elif conseq == 'stop retained':
                print(f'Variant at position {position}: Stop retained. No change to amino acid sequence. Skipped', flush= True)
            else:
                print(f'Variant at position {position}: Unknown consequenceType. Skipped.', flush= True)
                continue
            
            # add to dictionary
            fasta_seqs.append(SeqRecord(Seq(''.join(seq_new)), id= f've|{accession}_{desc}|{name}', description= f'{protname} OS={organism} OX={taxid} GN={genename} PE={protein_existence} SV={sequence_version}'))

        # write variant sequences to fasta
        with open(f'{outdir}/{accession}_{protname}_variants.fasta', 'w') as output:
            SeqIO.write(fasta_seqs, output, 'fasta')

        print(f'>>> Fasta file of variants for {os.path.basename(filepath)} successfully created.', flush= True)

    print('################## Completed! ##################', flush= True)

    return

###########################################################
###---------------------- TKINTER ----------------------###
###########################################################

import tkinter as tk
import os
from tkinter import filedialog, Listbox, Label, Entry, Button, messagebox, Text

class ConsoleRedirector:

    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, message):
        self.text_widget.insert(tk.END, message)
        self.text_widget.see(tk.END)

    def flush(self):
        pass

def browse_input_files():
    global file_list
    file_paths = filedialog.askopenfilename(title= "Select Input XML File", filetypes= [("XML files", "*.xml")], multiple= True)
    if file_paths:
        file_list = []
        file_listbox.delete(0, tk.END)
        for file in file_paths:
            if file not in file_list: # to avoid duplicates
                file_list.append(file)
                file_listbox.insert(tk.END, file)

def browse_output_folder():
    folder_path = filedialog.askdirectory(title= "Select Output Folder", initialdir= '.', mustexist= True)
    if folder_path:
        output_entry.config(state= "normal")
        output_entry.delete(0, tk.END)
        output_entry.insert(0, folder_path)
        output_entry.config(state= "readonly") 

def submit():
    output_folder = output_entry.get()
    if not file_list or not output_folder:
        messagebox.showerror("Error", "Both input file and output folder must be set.")
        return
    # Process function
    func_thread = threading.Thread(target= fasta_from_variantsXML, args= (file_list, output_folder), daemon= True)
    func_thread.start()

### main window
root = tk.Tk()
root.title("XML2FASTA")
#root.geometry("600x200")  # window size

file_list = []  # store selected file paths
output_folder = tk.StringVar() # store output folder path 

### select input file
Label(root, text="Input XML File(s):").pack()
btn_select_files = Button(root, text= "Browse Input File(s)", command= browse_input_files)
btn_select_files.pack(pady= 5)

file_listbox = Listbox(root, width=80, height=10)
file_listbox.pack(padx= 10, pady= 5)

### select output folder
Label(root, text="Output Folder:").pack()
btn_select_output = Button(root, text= "Browse Output Folder", command= browse_output_folder)
btn_select_output.pack(pady=5)

# **Output Folder Display**
output_entry = Entry(root, textvariable= output_folder, width= 80) #state= "readonly")
output_entry.pack(padx= 10, pady= 5)

# set initial folder (here)
output_entry.insert(0, os.path.abspath(os.getcwd())) # current working directory
output_entry.config(state= "readonly")

### Console Output Box
Label(root, text= "Console Output:").pack()
console_output = Text(root, height= 10, width= 80, state= "normal")
console_output.pack(padx= 10, pady= 5)

# redirect prints
sys.stdout = ConsoleRedirector(console_output)

### submit button
submit_button = Button(root, text= "Submit", command= submit)
submit_button.pack(pady= 20)

### Close Button
close_button = Button(root, text= "Close", command= root.destroy)
close_button.pack(pady= 10)

### run script
root.mainloop()