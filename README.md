# ConformationalChange
Matlab programs to detect local and large-scale conformational changes in RNA 3D structures

<h3>Dependencies</h3>

The ConformationalChange programs require a current release of FR3D, available at <a href="https://github.com/BGSU-RNA/FR3D/tree/dev">Github FR3D repository</a>.
If you want to use R3D Align to align 3D structures which have less than near-complete sequence identity, you will also need release 1.0 of R3D align, available at <a href="https://github.com/BGSU-RNA/R3DAlign/releases/tag/v1.0">Github R3D Align repository</a>.

We suggest placing the FR3D, R3DAlign, and ConformationalChange program files in three separate folders inside of a Github folder, and not to edit these programs directly so that your edits are not overwritten by downloading a new version of the repositories.

We suggest making a separate FR3D folder to be the working directory for Matlab and Octave.
Within this folder, you can make two sub-folders, one called PDBFiles, to store downloaded RNA 3D structure files, and the other called PrecomputedData, to store binary files created by FR3D to store its annotations and speed up processing the next time a 3D structure file is needed.

<h3>Setting paths in Matlab</h3>
It is useful to set the working directory to the FR3D folder, since FR3D will look in specific places for 3D structure files and pre-annotated versions of 3D structure files.  To set the working directory in Matlab, you can use the dropdown menu in the menu bar or use a command such as:

        cd('C:\Users\username\Documents\FR3D\');

You may want to configure your system so that Matlab starts in the desired working directory.

To tell Matlab where to find the ConformationalChange, FR3D, and R3D Align program files, choose File, Set Path, and add the following directories one by one, then click Save.

        C:\Users\username\Documents\GitHub\ConformationalChange
        C:\Users\username\Documents\GitHub\FR3D\FR3DSource
        C:\Users\username\Documents\GitHub\R3DAlign\R3DAlign
        C:\Users\username\Documents\FR3D\
        C:\Users\username\Documents\FR3D\PrecomputedData
        C:\Users\username\Documents\FR3D\PDBFiles


<h3>Setting paths in Octave</h3>
It is useful to set the working directory to the FR3D folder, since FR3D will look in specific places for 3D structure files and pre-annotated versions of 3D structure files.  To set the working directory in Octave, use a command such as:

        cd('C:\Users\username\Documents\FR3D\');

You may want to configure your system so that Octave starts in the desired working directory.

Next, tell Octave where to find the ConformationalChange, FR3D, and R3D Align program files by adding their directories to the path, using commands such as:

        addpath('C:\Users\username\Documents\GitHub\ConformationalChange');
        addpath('C:\Users\username\Documents\GitHub\FR3D\FR3DSource');
        addpath('C:\Users\username\Documents\GitHub\R3DAlign\R3DAlign');
        addpath('C:\Users\username\Documents\FR3D\');
        addpath('C:\Users\username\Documents\FR3D\PrecomputedData');
        addpath('C:\Users\username\Documents\FR3D\PDBFiles');

You might want to save all of these commands in a text file called SetPath.m, store it in the FR3D folder, and run SetPaths at the Octave prompt when you start Octave.

Note that Octave runs more slowly than Matlab.  Be patient.  Also, try as we might, we cannot suppress all Octave warning messages.

<h3>Running the examples</h3>

In Matlab and Octave, run programs by typing the names of the files at the command prompt, without the .m extension.

<h3>Modifying the examples or analyzing different files</h3>

We suggest copying the Example .m files from the ConformationalChange repository into the FR3D folder and editing the files there.
Each Example file is well documented to show how to change options, and it should be straightforward to make changes for different 3D structure files.
One good place to browse available RNA