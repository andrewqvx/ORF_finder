#!/bin/bash 
rm *.txt                                                                # removes text files in directory before starting program
echo input \(GID\)\:                                                    # prints "input GID"
read GID                                                                # stores user input value as a variable under the name "GID"
#GID=100226.110
#obtain file 
#wget ftp://anonymous:password@ftp.bvbrc.org/genomes/$GID/$GID.fna      # obtains fna file from bv-brc database

echo creating forward strand...
#turn fna file into one line fasta sequence 
cat $GID.fna | tr '\n' '#' | cut -d '#' -f 2- | tr -d '#' > for_strand.txt # replaces a '#' at the end of every new line | uses # as a custom delimiter and keeps everything in the second                                                                              field and beyond | deletes #; results in a one line fasta sequence

echo creating reverse complement strand...
rev for_strand.txt | tr "atgc" "tacg" > rev_strand.txt                  # reverse command reverses all characters in the text file, regardless if they are on a new line. IMPORTANT to                                                                               utilize this command and not "tac" as this will result in reversing the order of lines where as "rev" will                                                                                 reverse the order of characters | "tr" command takes the following characters in the first arguement and replaces                                                                          them with characters in the second argument
echo creating open reading frames...
#create open reading frames
sed 's/.\{3\}/&\_/g' for_strand.txt > for1.txt                          # sed command uses the substitiution operation "s/"
cut -c 2- for_strand.txt | sed 's/.\{3\}/&\_/g' > for2.txt              # '.' = a wildcard that searches for any three characters
cut -c 3- for_strand.txt | sed 's/.\{3\}/&\_/g' > for3.txt              # \{3\} = reg. expression instructing every three times
sed 's/.\{3\}/&\_/g' rev_strand.txt > rev1.txt                          # the 2 lines above match any three characters
cut -c 2- rev_strand.txt | sed 's/.\{3\}/&\_/g' > rev2.txt              # & matches the first argument of the sed command
cut -c 3- rev_strand.txt | sed 's/.\{3\}/&\_/g' > rev3.txt              # '_' is appended after the first argument
                                                                        # the 2 lines above take the first argument and replace the same found argument with an underscore after
                                                                        # /g is a flag that = global. Will perform this command for every occuence within a line...not just the first occurence
echo finding open reading frames...

#find open reading frames in each reading frame

sed -E 's/taa|tag|tga/&\n/g' for1.txt | sed 's/atg/*&/' | grep '*' | cut -d '*' -f 2 | grep -E 'taa|tag|tga' | tr -d '_' | grep -E '.{75,}' | sed 's/^atg/*&/' | awk 'BEGIN { orf=0 } {gsub(/\*/, ">Fr1." (++orf) "\n")} 1' > orf.txt
sed -E 's/taa|tag|tga/&\n/g' for2.txt | sed 's/atg/*&/' | grep '*' | cut -d '*' -f 2 | grep -E 'taa|tag|tga' | tr -d '_' | grep -E '.{75,}' | sed 's/^atg/*&/' | awk 'BEGIN { orf=0 } {gsub(/\*/, ">Fr2." (++orf) "\n")} 1' >> orf.txt
sed -E 's/taa|tag|tga/&\n/g' for3.txt | sed 's/atg/*&/' | grep '*' | cut -d '*' -f 2 | grep -E 'taa|tag|tga' | tr -d '_' | grep -E '.{75,}' | sed 's/^atg/*&/' | awk 'BEGIN { orf=0 } {gsub(/\*/, ">Fr3." (++orf) "\n")} 1' >> orf.txt
sed -E 's/taa|tag|tga/&\n/g' rev1.txt | sed 's/atg/*&/' | grep '*' | cut -d '*' -f 2 | grep -E 'taa|tag|tga' | tr -d '_' | grep -E '.{75,}' | sed 's/^atg/*&/' | awk 'BEGIN { orf=0 } {gsub(/\*/, ">Fr-1." (++orf) "\n")} 1' >> orf.txt
sed -E 's/taa|tag|tga/&\n/g' rev2.txt | sed 's/atg/*&/' | grep '*' | cut -d '*' -f 2 | grep -E 'taa|tag|tga' | tr -d '_' | grep -E '.{75,}' | sed 's/^atg/*&/' | awk 'BEGIN { orf=0 } {gsub(/\*/, ">Fr-2." (++orf) "\n")} 1' >> orf.txt
sed -E 's/taa|tag|tga/&\n/g' rev3.txt | sed 's/atg/*&/' | grep '*' | cut -d '*' -f 2 | grep -E 'taa|tag|tga' | tr -d '_' | grep -E '.{75,}' | sed 's/^atg/*&/' | awk 'BEGIN { orf=0 } {gsub(/\*/, ">Fr-3." (++orf) "\n")} 1' >> orf.txt

# sed command uses a substitution operation to search globally in desired txt file for any of the stop codons and replace the stop codon with the same text followed by a new line. E option allows the use of special character "|" expressing an 'or' arguement
# sed command uses a substitution to search for a start codon and place a '*' in front of the first atg found on each line
# grep command prints the lines that only contain the '*' within the line, eliminating any lines do not contain atg
# cut command uses '*' as custom delimiter and keeps fields that come after the first '*'
# grep command keeps lines that contain one of the three stop codons...eliminates the lines that do not meet the criteria for an orf
# tr command deletes all the underscores within the txt file
# grep searches for any lines with 75+ characters. -E is an option to allow the use of metacharacters... {} in this command. '.' = any charcter EXCEPT a new line character.
# sed command uses the substiitution operation to search the any line that begins with 'atg' and replace the pattern with '*' followed by the same search pattern. This occurs once on every new line.
#awk command contains an block before awk command begins processing the argument(input).
    # {orf=0} is a variable that will count each line that meets an argument)
    # gsub = function called gsub (global substitution) that searches for a pattern and replaces it
    # \* (++counter) "\n") = searches for '*'
    #  >Fr-3. = replacement pattern; replaces every '*' with this pattern.
    # (++counter) "\n" = orf counter will print count number directly after replacement pattern and then add a newline after.
    # 1 prints the modified line
    
    #!!awk commmand: will process input lines and replace each * with a sequence like >Frame1.1, >Frame1.2, and so on, incrementing the value for each occurrence of *. The modified lines are then printed as output.
echo creating proteins from ORFs...

####### loop time

IFS=$'\n'                                                               # Internal Field Separator: every new line is a new field that is to be processed
for ORF in $(cat orf.txt)                                               # for loop that takes each line and stores it in an arbitrary variable named ORF
do

length=${#ORF}                                                          # calculates length of each ORF iteration sequence and stores it as a variable named 'length'

    if [[ $ORF == *">"* ]]                                              # if statement that appends metadata line to a temp file, ultimately creating the structure for the protein fasta file, conditional statement checks if the value of the variable $ORF contains the > character, and it evaluates to true if it does
    then
    echo $ORF >> temp2.txt
    else
    : #null command - no operation
    fi #ends if block (backwards if); to not process if statement further
    
    for (( n=0; n < length; n+=3 ))                                     # nested loop! lines that do NOT contain metadata get passed through this counter loop. Takes each ORF iteration and                                                                     increments three character until 'n' reaches the end of the ORF length. Each increment of three characters will be                                                                     further processed.
    do
    
    codon="${ORF:n:3}"                                                  # extracts a sequence from ORF from the nth position and cut to a length of three characters starting at position n.                                                                     The extracted sequence is stored in the codon variable which is further processed by 'if' statements
    
        if [[ $ORF == *">"* ]]                                          # this 'if' statment instructs the for loop to pass or ignore the metadata lines and continue to the next loop
        then                                                            # iteration. REMEMBER that we already appended this metadata line to a temp text file in the first 'for' loop. We
        continue                                                        # need this if statment in this nested loop so the counter loop doesn't try to process  metadata line
        fi
        
        if [ "${codon}" == "atg" ]                                      # if codon variable equals the codon's listed in the respective 'if' statement, then it will print the respective
        then                                                            # amino acid letter
        echo "M" >> temp.txt
        
        elif [ "${codon}" == "ttt" ] || [ "${codon}" == "ttc" ]         #elif used when if statement is false or doesn't match with previous if/elif statement
        then
        echo "F" >> temp.txt

        elif [ "${codon}" == "tta" ] || [ "${codon}" == "ttg" ] || [ "${codon}" == "ctt" ] || [ ""${codon}"" == "ctc" ] || [ "${codon}" == "cta" ] || [ ""${codon}"" == "ctg" ]
        then
        echo "L" >> temp.txt
        
        elif [ "${codon}" == "att" ] || [ "${codon}" == "atc" ] || [ "${codon}" == "ata" ]
        then
        echo "I" >> temp.txt
        
        elif [ "${codon}" == "gtt" ] || [ "${codon}" == "gtc" ] || [ "${codon}" == "gta" ] || [ "${codon}" == "gtg" ]
        then
        echo "V" >> temp.txt
        
        elif [ "${codon}" == "tct" ] || [ "${codon}" == "tcc" ] || [ "${codon}" == "tca" ] || [ "${codon}" == "tcg" ] || [ "${codon}" == "agt" ] || [ "${codon}" == "agc" ]
        then
        echo "S" >> temp.txt
        
        elif [ "${codon}" == "cct" ] || [ "${codon}" == "ccc" ] || [ "${codon}" == "cca" ] || [ "${codon}" == "ccg" ]
        then
        echo "P" >> temp.txt
        
        elif [ "${codon}" == "act" ] || [ "${codon}" == "acc" ] || [ "${codon}" == "aca" ] || [ "${codon}" == "acg" ]
        then
        echo "T" >> temp.txt
       
        elif [ "${codon}" == "gct" ] || [ "${codon}" == "gcc" ] || [ "${codon}" == "gca" ] || [ "${codon}" == "gcg" ]
        then
        echo "A" >> temp.txt
        
        elif [ "${codon}" == "tat" ] || [ "${codon}" == "tac" ]
        then
        echo "Y" >> temp.txt
        
        elif [ "${codon}" == "cat" ] || [ "${codon}" == "cac" ]
        then
        echo "H" >> temp.txt
        
        elif [ "${codon}" == "caa" ] || [ "${codon}" == "cag" ]
        then
        echo "Q" >> temp.txt
        
        elif [ "${codon}" == "aat" ] || [ "${codon}" == "aac" ]
        then
        echo "N" >> temp.txt
        
        elif [ "${codon}" == "aaa" ] || [ "${codon}" == "aag" ]
        then
        echo "K" >> temp.txt
        
        elif [ "${codon}" == "gat" ] || [ "${codon}" == "gac" ]
        then
        echo "D" >> temp.txt
        
        elif [ "${codon}" == "gaa" ] || [ "${codon}" == "gag" ]
        then
        echo "E" >> temp.txt
        
        elif [ "${codon}" == "tgt" ] || [ "${codon}" == "tgc" ]
        then
        echo "C" >> temp.txt
        
        elif [ "${codon}" == "tgg" ]
        then
        echo "W" >> temp.txt
        
        elif [ "${codon}" == "cgt" ] || [ "${codon}" == "cgc" ] || [ "${codon}" == "cga" ] || [ ""${codon}"" == "cgg" ] || [ "${codon}" == "aga" ] || [ ""${codon}"" == "agg" ]
        then
        echo "R" >> temp.txt
        
        elif [ "${codon}" == "ggt" ] || [ "${codon}" == "ggc" ] || [ "${codon}" == "gga" ] || [ "${codon}" == "ggg" ]
        then
        echo "G" >> temp.txt
        
        elif [ "${codon}" == "taa" ] || [ "${codon}" == "tag" ] || [ "${codon}" == "tga" ]
        then
        echo "*" >> temp.txt
        
        else
        echo "error" >> temp.txt
        
        fi
        
    tr -d [:space:] < temp.txt >> temp2.txt                                         # deletes all blank space essentially creating one line from single characte lines that the for loop                                                                                     created
    
    > temp.txt                                                                      # removes text from temp file, before counter loop begins again.
    
    done

                                              
done
sed 's/\*/&\n/g' temp2.txt >> amino.txt                                             # sed command uses substitution operation to search for '*' globally and replace it with the same                                                                                            character followed by a new line. In essensce, this will the structure for protein fasta file.

for1_num=$(grep '>Fr1' orf.txt | wc -l)                                             # create variables for each reading frame that will obtain orf count
for2_num=$(grep '>Fr2' orf.txt | wc -l)
for3_num=$(grep '>Fr3' orf.txt | wc -l)
rev1_num=$(grep '>Fr-1' orf.txt | wc -l)
rev2_num=$(grep '>Fr-2' orf.txt | wc -l)
rev3_num=$(grep '>Fr-3' orf.txt | wc -l)


echo number of open reading frames \in Fr\1\: $for1_num                             # print number of orfs found in each of the 6 reading frames
echo number of open reading frames \in Fr\2\: $for2_num
echo number of open reading frames \in Fr\3\: $for3_num
echo number of open reading frames \in Fr\-1\: $rev1_num
echo number of open reading frames \in Fr\-2\: $rev2_num
echo number of open reading frames \in Fr\-3\: $rev3_num

orf_num=$(grep '>' orf.txt | wc -l)                                                 # print total amount of open reading frames
echo total number of open reading frames found\: $orf_num


protein_num=$(grep '>' amino.txt | wc -l)                                           # creating variable that counts the total number of proteins created from ORFs
echo total number of proteins found\:$protein_num

awk '/^>/ {if (seq) {print seq; seq=""; print ""} print; next} {seq = seq $0} END {if (seq) print seq}' orf.txt | fold -w 60 > orf.fasta
awk '/^>/ {if (seq) {print seq; seq=""; print ""} print; next} {seq = seq $0} END {if (seq) print seq}' amino.txt | fold -w 60 > protein.fasta


# '/^>/ {if (seq) {print seq; seq=""; print ""} print; next}:
# /^>/ = awk pattern that matches lines starting with >
    # if (seq) {print seq; seq=""; print ""} = if the seq variable is not empty, it prints its contents (sequence data), sets seq to an empty string (resetting it), and then prints an empty  line creating a space between sequences and metadata headers
    # print = prints the header line starting with >.
    # next = skips the rest of the awk processing for this line and moves to the next line in the input
    
# {seq = seq $0} = action block for lines that dont start with > appends the current line to the variable (sequence in this instance

# END {if (seq) print seq} = block is executed at the end of processing, after all lines have been processed
    # if (seq) print seq = if seq variable is not empty, it means theres unprinted sequence data remaining, so it prints it
    
#!!awk command reads a FASTA-formatted file, separates sequence headers from sequence data, and prints the sequences with each sequence on a separate line. It handles the case where sequences are split across multiple lines. When it encounters a new sequence header (line starting with >), it prints the previous sequence and starts accumulating the new sequence. It prints a blank line between sequences for clarity.

