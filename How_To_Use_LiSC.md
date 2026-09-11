# Connecting to LiSC — Workshop Setup Guide

This course uses the **Life Science Compute Cluster (LiSC)** of the University of Vienna. A LiSC account has already been created for you as part of the course group `2026s_evolvienna`. You were **not** sent a welcome email or a temporary password directly — instead, follow the steps below to activate your account before the course starts.

## 1. Allow your network to reach LiSC

Before anything else will work, go to **https://lisc.univie.ac.at/firewall/** and enter your LiSC username to allow your current network connection through the firewall.

This access is always **temporary** whenever you are outside the University of Vienna network — it lasts about **12 hours**, so you'll need to repeat this step again after that if you're still off-network (e.g. on home Wi-Fi, a hotel network, or a VPN).

**During the workshop, you will not need to do this at all if you connect via eduroam at the department** — eduroam there is already inside the trusted network range.

## 2. Your LiSC username

Your LiSC username was assigned by the cluster admin and shared with you by the course organizers via email (usually your last name, e.g. `smith`).

## 3. Activate your account and set a password

1. Go to **https://account.lisc.univie.ac.at**
2. Click **"Forgot password?"**
3. Enter your LiSC username and the email address you used to **register for this workshop**.
4. You will receive an email with a password-reset link — follow it to set your password.

Notes:
- Your workshop account expires automatically on **2026-09-30**.

## 4. Log in

Connect via SSH to one of the two login nodes:

```bash
ssh <your-username>@login01.lisc.univie.ac.at
# or
ssh <your-username>@login02.lisc.univie.ac.at
```

- **Mac / Linux:** use the built-in Terminal app.
- **Windows:** use the Windows Terminal / PowerShell (OpenSSH is built in), or a tool like PuTTY or MobaXterm.

## 5. Where to work during the course

- **Shared course folder** (read/write for everyone in the course):
  `/lisc/data/scratch/course/2026s_evolvienna`
  This folder is **not backed up** and will be **deleted automatically on 2026-09-30**. Please copy anything you want to keep elsewhere before then.
- Your personal home directory is also available for your own copies and exercises.

### Create your own subfolder — don't work directly in the shared folder

Everyone in the course has write access to the shared folder above, so please don't leave loose files there. Instead, on your first login, create your own subfolder named after your LiSC username and lock it so nobody else can write into it:

```bash
cd /lisc/data/scratch/course/2026s_evolvienna
mkdir <your-username>
chmod 700 <your-username>
```

`chmod 700` makes the folder fully private: only you (the owner) can read, write, or enter it — nobody else in the course, not even other participants or TAs, can access it. Do all your work for the course inside this folder.

## 6. Quick orientation once you're logged in

```bash
pwd                 # show where you are
ls                  # list files in the current folder
cd some_folder      # move into a folder
module avail        # list available software
module load R       # load a piece of software (e.g. R) into your session
exit                # log out
```

## 7. Software you need for this course

None of the software below is loaded by default when you log in — you load it into your session with `module load`. Please check this **before** the workshop starts:

```bash
git clone https://github.com/popgenomics/workshop_vienna.git
cd workshop_vienna
module load R-abcrf/1.9 IQ-TREE/3.1.3 ASTER/1.25-foss-2025a SciPy-bundle/2023.07-gfbf-2023a
bash server/check_prerequisites.sh
```

If the last line printed is `Blocking errors: 0`, you're ready for the workshop. If not, contact us before the course starts so we can help you fix it.

| Purpose | Module to load | Provides |
|---|---|---|
| R + required packages (tidyverse, abcrf) | `R-abcrf/1.9` | R, tidyverse, abcrf 1.9 |
| Phylogenetics (Day 2) | `IQ-TREE/3.1.3` | `iqtree3` |
| Species tree (Day 2) | `ASTER/1.25-foss-2025a` | `astral` |
| Python + pandas (Day 2 post-processing) | `SciPy-bundle/2023.07-gfbf-2023a` | Python 3.11, pandas 2.0.3 |

**Important:** use `module load R-abcrf/1.9`, not the generic `module load R` — the plain `R` module has a different, incompatible version of the `abcrf` package. You'll need to run these `module load` commands again every time you start a new login session.

Compiling the Aphid tool (Day 2) uses the system's default `gcc` — no module needed for that.

## 8. Transferring files

Use SCP, SFTP, or rsync to the same login nodes.

### Uploading files to LiSC

```bash
scp /path/to/my_file.txt <your-username>@login01.lisc.univie.ac.at:~/

rsync -avP /path/to/my_file.txt <your-username>@login01.lisc.univie.ac.at:~/
```

### Downloading files from LiSC

```bash
scp <your-username>@login01.lisc.univie.ac.at:/path/to/remote_file.txt /path/to/local_folder/

rsync -avP <your-username>@login01.lisc.univie.ac.at:/path/to/remote_file.txt /path/to/local_folder/
```

### Using an IDE instead

You don't have to use the command line for file transfer. An IDE such as **VS Code** with the **Remote - SSH** extension lets you connect directly to a LiSC login node, browse the remote filesystem in a normal file explorer sidebar, and drag-and-drop files to upload or download.

## 9. Getting help

- Documentation: https://wiki.lisc.univie.ac.at (see "Getting started" and "Howto")
- Guided tutorial (assemble/annotate a virus genome): https://wiki.lisc.univie.ac.at/access/gettingstarted/tutorial
- Monthly intro webinar (LiSC structure, SLURM, Unix/bash basics), second Thursday of the month, 15:00-16:30: https://wiki.lisc.univie.ac.at/support/webinar
- Technical problems / bugs: https://lisc.univie.ac.at/helpdesk

## Terms of use

Your LiSC account is subject to the LiSC terms of use (linked from https://lisc.univie.ac.at). Keep your password private and never share your account with anyone else.
