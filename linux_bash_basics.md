# Introduction to Linux / Bash for Bioinformatics

*A quick-start guide for biologists working on a Linux computing cluster*

---

## 1. Why Linux and the Command Line?

Most bioinformatics tools (BLAST, samtools, GATK, BWA, etc.) and computing clusters run on **Linux**. Instead of clicking icons, you interact with the computer by typing commands in a **terminal** — this is called the **command line** or **shell** (here, we use **Bash**, the most common shell).

It looks intimidating at first, but with about 15 commands you can already do 90% of what you need.

---

## 2. Opening a Terminal and Connecting to the Cluster

To connect to a remote cluster, you typically use `ssh` (Secure Shell):

```bash
ssh your_username@cluster_address
```

Example:
```bash
ssh jdupont@cluster.university.fr
```

You will be asked for your password (note: when typing, **no characters appear on screen** — this is normal, just type and press Enter).

---

## 3. Understanding the Prompt

Once connected, you'll see something like:

```bash
jdupont@cluster:~$
```

This tells you:
- `jdupont` → your username
- `cluster` → the machine name
- `~` → your current location (here, your **home directory**)
- `$` → waiting for your command

---

## 4. Navigating the File System

Linux organizes files in a **tree structure**, starting from the root `/`.

| Command | What it does |
|---|---|
| `pwd` | **P**rint **W**orking **D**irectory — shows where you are |
| `ls` | **L**i**s**t files/folders in the current directory |
| `ls -l` | List with details (size, date, permissions) |
| `ls -lh` | Same, but with **h**uman-readable sizes (KB, MB, GB) |
| `ls -a` | List including hidden files (starting with `.`) |
| `cd folder_name` | **C**hange **D**irectory — move into a folder |
| `cd ..` | Move up one level (to the parent folder) |
| `cd ~` or `cd` | Go back to your home directory |
| `cd /path/to/somewhere` | Go directly to a specific location (absolute path) |

### Example session
```bash
pwd
# /home/jdupont

ls -lh
# drwxr-xr-x  2 jdupont group  4.0K Jan 10 09:15 raw_data
# -rw-r--r--  1 jdupont group  2.3M Jan 10 09:20 genome.fasta

cd raw_data
pwd
# /home/jdupont/raw_data
```

> 💡 **Tip:** Press `Tab` while typing a file/folder name — Bash will try to auto-complete it for you. This saves time and avoids typos!

---

## 5. Creating, Copying, Moving, and Deleting

| Command | What it does |
|---|---|
| `mkdir new_folder` | Create a new folder |
| `touch file.txt` | Create an empty file |
| `cp file.txt copy.txt` | Copy a file |
| `cp -r folder1 folder2` | Copy a folder and its contents (recursive) |
| `mv file.txt new_name.txt` | Rename a file |
| `mv file.txt some_folder/` | Move a file into a folder |
| `rm file.txt` | Delete a file (⚠️ no trash bin, it's permanent!) |
| `rm -r folder_name` | Delete a folder and everything inside it (⚠️ **be very careful**) |

> ⚠️ **Warning:** There is no "Recycle Bin" in the terminal. `rm` deletes files immediately and permanently. Always double-check the file/folder name before pressing Enter, especially with `rm -r`.

---

## 6. Looking Inside Files

Bioinformatics files (FASTA, FASTQ, VCF, CSV...) are often plain text, so you can preview them without opening a heavy editor.

| Command | What it does |
|---|---|
| `cat file.txt` | Print the **whole** file to the screen |
| `head file.txt` | Show the **first 10 lines** |
| `head -n 20 file.txt` | Show the first 20 lines |
| `tail file.txt` | Show the **last 10 lines** |
| `less file.txt` | Open the file for scrolling (press `q` to quit) |
| `wc -l file.txt` | Count the number of lines |

### Example: peeking at a FASTA file
```bash
head -n 4 genome.fasta
# >chr1
# ATGCGTACGTTAGCATGC...
# >chr2
# GGCATTAGCATGCATGCA...
```

---

## 7. Searching Inside Files

| Command | What it does |
|---|---|
| `grep "pattern" file.txt` | Find lines containing "pattern" |
| `grep -c "pattern" file.txt` | Count matching lines |
| `grep -v "pattern" file.txt` | Show lines that **do NOT** match |
| `grep -i "pattern" file.txt` | Case-insensitive search |

### Example: counting sequences in a FASTA file
```bash
grep -c ">" genome.fasta
# 24   (number of sequences, since each starts with ">")
```

---

## 8. Combining Commands: the Pipe `|`

The pipe `|` sends the output of one command as input to the next. This lets you build powerful one-liners.

```bash
grep ">" genome.fasta | wc -l
# counts how many sequences are in the FASTA file
```

```bash
ls -lh *.fastq | head -5
# shows details of the first 5 fastq files
```

---

## 9. Wildcards: Selecting Multiple Files

| Symbol | Meaning |
|---|---|
| `*` | Any number of characters |
| `?` | Exactly one character |

```bash
ls *.fastq        # all files ending in .fastq
ls sample_?.txt   # sample_1.txt, sample_2.txt, etc. (not sample_10.txt)
```

---

## 10. Redirecting Output: `>` and `>>`

| Symbol | Meaning |
|---|---|
| `>` | Save output to a file (overwrites it!) |
| `>>` | Append output to the end of a file |

```bash
grep ">" genome.fasta > sequence_ids.txt
echo "Analysis done" >> logfile.txt
```

---

## 11. File Permissions (Quick Overview)

```bash
ls -l script.sh
# -rwxr-xr-x  1 jdupont group  512 Jan 10 script.sh
```

The `rwx` triplets mean **r**ead, **w**rite, e**x**ecute — for (1) the owner, (2) the group, (3) everyone else.

To make a script executable:
```bash
chmod +x script.sh
./script.sh
```

---

## 12. Running Jobs on a Cluster (Just the Idea)

Clusters use a **job scheduler** (e.g., Slurm, PBS, SGE) so many users can share resources. You typically don't run heavy analyses directly — you submit a **job**:

```bash
sbatch my_analysis.sh      # Slurm example
squeue -u jdupont           # check the status of your jobs
```

*(The exact commands depend on your specific cluster — this will likely be covered in a dedicated session.)*

---

## 13. Getting Help

| Command | What it does |
|---|---|
| `man command_name` | Manual page for a command (press `q` to quit) |
| `command_name --help` | Short help message |

```bash
man ls
grep --help
```

---

## 14. Cheat Sheet — The Essentials

```bash
pwd                 # where am I?
ls -lh               # what's here?
cd folder/           # go into folder
cd ..                # go back up
mkdir results        # create a folder
cp file.txt dest/    # copy a file
mv old.txt new.txt   # rename/move
rm file.txt          # delete (careful!)
head -n 5 file.txt   # peek at a file
grep "word" file.txt # search in a file
command1 | command2  # chain commands
command > out.txt    # save output
Tab                  # auto-complete
Ctrl + C             # cancel a running command
Ctrl + R             # search command history
```

---

## 15. Practice Exercise

1. Connect to the cluster with `ssh`.
2. Create a folder called `workshop_test`.
3. Move into it with `cd`.
4. Create an empty file called `notes.txt`.
5. Write "Hello Linux" into it using `echo "Hello Linux" > notes.txt`.
6. Display its content with `cat notes.txt`.
7. Copy it to `notes_backup.txt`.
8. List all files in the folder with `ls -l`.
9. Delete `notes_backup.txt`.

If you can do all of this, you're ready for the rest of the workshop! 🎉
