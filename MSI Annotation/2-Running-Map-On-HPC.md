# Running MAP
After following the previous instructions, you should have an input .txt file ready for input into the [MAP](https://sourceforge.net/p/mapmsi/wiki/MAP/) software! MAP was made as an R package. However, when I was testing it out for myself in R studio, one of the dependencies ("GSVA", or gene set variation analysis) that MAP is based on was too recent. In other words, MAP was created in 2020 and used the old version of the dependency that R studio no longer supports. To install the old version of the dependency, I would need to download the old version of R studio, which would be too much of a hassle.

To work around this problem, we will use our High-Performance Computing (HPC) system (aka supercomputer)! The primary advantage of the supercomputer is that it accomodates intensive computations with large files that are not possible on our own computers. However, the secondary advantage is that it provides a shared project space with large storage, where I can deposit files and packages into the project space for everyone to access WITHOUT having to download these for yourself!

In the HPC computing world, there is something called a container. Essentially, this is a sharable file/"package" of software containing all necessary libraries, environments, and dependencies of a project. The most commonly used container system is known as Docker, but I will be using Singularity container system because it was created specifically for HPC systems. To learn more about Singularity and containers, here is a [source](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) you can explore. I also found some additional [youtube](https://www.youtube.com/watch?v=_KhbwXqk0Bk) [videos](https://www.youtube.com/watch?v=xty42A05Wg0) that may be helpful (this is 100% optional and not needed for the instructions of this project).

A Singularity container is stored in the form of a .sif file. In the container I have created, I was able to install MAP, the old version of R, and thus the old version of the GSVA dependency. Side note: one major advantage of containers is that you may store versions of software that is completely different from the rest of your computer's environment.

The following instructions will teach you how to access and use this container to run MAP!

## Table of Contents
1: Supercomputer setup and information to know
* 1.1: Setting up and accessing your account
* 1.2: Directories and paths to know
* 1.3: Files needed and where to store them

2: Setting up R script

3: Setting up batch script

4: Running batch script

5: Accessing MAP output

## Step 1: Supercomputer setup and information to know
First, you will need one of the lab leadership members to create an Expanse supercomputer account for you (this should have been done when you first joined, but if not let someone know!). 

### 1.1: Setting up and accessing your account
For detailed instructions to access and log in to your supercomputer account, please work through the instructions I have compiled in my lab notebook, under [Supercomputer- 1: Access and Login](https://github.com/alfalfacow/Lab-Notebook/blob/main/Supercomputer/1_Access_and_Login.md). 

### 1.2: Directories and paths to know
There are three main "paths" within the supercomputer account that you will be working with. Paths tell you about folder hierarchy and where to navigate (small example: /home/downloads/hi.txt tells you that hi.txt is a file in the downloads folder within the home directory).

The first is your "home" directory, which looks like /home/user/etc/etc (replace user with your own username). This has a limited storage and is where you may store more permanent things like environments and small important downloads.

Another main location to know is the "[SDSC Expanse Lustre file system](https://www.sdsc.edu/systems/expanse/user_guide.html)". The following two paths/directories are within this file system. 

The second "path" is your "lustre scratch" directory, which looks like /expanse/lustre/scratch/user/etc/etc (replace user with your own username). Unlike the home directory, this is meant to be temporary storage for your projects (inactive files are apparently purged/deleted 90 days after our project expires). It also holds LOT more storage (I believe up to 1 terabyte!!). Other lab members in the same project do not have access to this directory.

The third "path" is the "lustre project" directory, which looks like /expanse/lustre/scratch/project/user/etc/etc (replace user with your own username, and project with csd670 (our current project ID)). This is a shared space for everyone in the project. Everyone has their own folder (under "user" in the path) that they can add files to, and all other lab members can access items in this folder.

This last folder (shared project folder) is where I have stored my Singularity container .sif file for you to use. For your reference, the path to this file is:
```
/expanse/lustre/projects/csd670/akao1/MAP/map_msi.sif
```
### 1.3: Files needed and where to store them
To run MAP, you need to prepare 4 files:
* .sif container file (see above)
* Input .txt file (log2 transformed gene expression matrix) (section 1 of these instructions, store these within your supercomputer account, whether that be in your home directory or your lustre scratch directory)
* R script to run MAP (see step 2 below)
* Batch script to submit the final command (see step 3 below)

For the purposes of this project, I encourage you to create and work within a single directory out of one of the 3 paths above. That way, all of your files and results will be centralized in a single folder. Examples are shown below:
```
#Using home directory
cd /home/user #changes current directory to this one, although this should be the default
mkdir /home/user/MAP #creates new folder at the designated path

#Using expanse lustre scratch
cd /expanse/lustre/scratch/user/temp_project
mkdir /expanse/lustre/scratch/user/temp_project/MAP

#Using expanse lustre project
cd /expanse/lustre/projects/csd670/user
mkdir/expanse/lustre/projects/csd670/user/MAP
```
Aside from the .sif container file, the remaining 3 files will be prepared by you and should be deposited in the new project directory of your choosing!

## Step 2: Setting up R script
Setting up the R script is very simple (only 4 lines!) You may do this within R studio (create new script and then export it as a .R file), or edit a google doc and export it as a .txt file (and then renaming the extension to .R). Personally I haven't tried the google doc way (not 100% sure if it's reliable) so I prefer doing this step within R studio. 

On the top right of R studio there should be green plus on top of a blank paper icon. Click "R script" (or just use the shift+command+N shortcut) to open a new script tab. Then, input the following code:
```
library(MAP)
testP = "/PATH/TO/YOUR/DESIRED/OUTPUT/DIRECTORY"
expData = "/PATH/TO/YOUR/MAP_input.txt"
runMAP(expData, testP)
```
testP is the folder/directory where you want to store your MAP results. In our case, we would store it in the project directory created in step 1.3.

expData is the path to the input MAP_input.txt file you created in the previous instructions!

As mentioned in Step 1, it is good practice to keep all work for a specific project in its own folder/directory. We can designate testP as a new output folder within the directory we chose before, and expData as the directory we chose before (where we stored the .txt input file).
```
#Example if we had used expanse lustre scratch as the project directory
testP = "/expanse/lustre/scratch/user/temp_project/MAP/output" #points to folder named "map_output" within the MAP project directory we created
expData = "/expanse/lustre/scratch/user/temp_project/MAP/MAP_input.txt" #points to specific file MAP_input.txt within the MAP project directory we created

#OR

#Example if we had used expanse lustre project as the project directory
testP = "/expanse/lustre/projects/csd670/user/MAP/map_output" #points to folder named "map_output" within the MAP project directory we created
expData = "/expanse/lustre/projects/csd670/user/MAP/MAP_input.txt" #points to specific file MAP_input.txt within the MAP project directory we created
```

Next, save your R script and rename it anything you want (I named it run_map.R). Finally, transfer this .R file into the project directory we created in step 1.3! You can do this easily through the SFTP tab on Termius (see ending section of my [notes](https://github.com/alfalfacow/Lab-Notebook/blob/main/Supercomputer/1_Access_and_Login.md)).

## 3: Setting up batch script
If this is your first time reading about what a "batch script" is, I highly recommend you go through my [notes on supercomputer in general where I explain batch jobs](https://github.com/alfalfacow/Lab-Notebook/blob/main/Supercomputer/3_Submitting_Batch_Jobs.md).

This is the template you will be editing for this specific batch script: [click this link](https://docs.google.com/document/d/1KexZKKBtnPk90qh5hP3uG_Mi06yx-gyQvXfIILS4pwQ/edit?usp=sharing). It is view only, so make a copy to edit. The contents have been copied below:

```
#!/bin/bash
#SBATCH --job-name=MAP_MSI
#SBATCH --output=MAP_output.txt
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL
#SBATCH -t 1:00:00
#SBATCH --mail-user=YOUR_EMAIL_HERE@ucsd.edu
#SBATCH --mail-type=all
#SBATCH -A CSD670
#SBATCH --mem=20G

echo "job start"

module load singularitypro

SIF=/expanse/lustre/projects/csd670/akao1/MAP/map_msi.sif
DATA_DIR=/expanse/lustre/projects/csd670/akao1/MAP
SCRIPT=$DATA_DIR/run_map.R

singularity exec --bind $DATA_DIR:$DATA_DIR $SIF Rscript $SCRIPT

echo "job done" 
```

There are several areas that you will need to edit yourself:
* #SBATCH --mail-user (set this to your own email)
* DATA_DIR= (set this to the project folder/directory that we set up in step 1.3)
* SCRIPT= (set this to $DATA_DIR/the_name_of_your_R_script.R)

After you are done editing these fields, save the google doc as a .txt file and transfer it to the project directory via Termius SFTP!

## 4: Running batch script
At this point, you should double check that all the files are in the correct place. Here is what the hierarchy should be like:

>Parent path/directory (example: /expanse/lustre.... etc etc )
->MAP (project directory name)
-->output (folder) #if this folder doesn't exist yet, the MAP script should create it so don't worry if it isn't there yet
-->MAP_input.txt (input file)
-->run_map.R (R script)
-->MAP_script.txt (batch script)

Once you have verified this, you are ready to submit the batch job!

First, we run a command to configure the google document format into the necessary format:
```
cd /path/to/your/project/directory #change current directory to your project directory (whichever one you made) to access the sbatch file
dos2unix MAP_script.txt
```
Then we can send off the script as a batch job using sbatch!
```
sbatch MAP_script.txt
```
You should recieve an email telling you when your job has started and stopped. When you navigate to the output folder in your project directory through the SFTP tab, you should find the output files (and thus the predicted MSI status of the sample)!!! 

Congratulations, you have successfully predicted the MSI status of a colorectal cancer tissue sample! :)
