# atlaslc

### Table of Contents:
* [Installation](#installation)
* [Usage](#usage)
* [Additional usage](#additional-usage)
* [Example commands](#example-commands)

### Description: 
Using a table of SN names, RA, and Dec, atlaslc logs into the ATLAS machines and retrieves the SN light curve data. You will need a table containing your **SN names and your username and password** for the ATLAS machines.

### Installation:
* Download all code.
* Create a directory that will house the source code you downloaded and the data files you'll download in the future.
* Within that directory, create and name two additional directories. One will be your **source directory**, and the other will be your **data directory**. 
* Move all code into your **source directory**.
* In the **source directory**, open `atlaslc.sourceme`.
* Add the following code to the `<atlaslc.sourceme>` file below the other `elif` statements: 
	* `elif [[ $HOSTNAME =~ ADDHOSTNAMEHERE ]]; then`
	* `export ATLASLC_SOURCEDIR=`
  * `export ATLASLC_DATA=`
  * `export PIPE_BATCH_SYSTEM=NONE`
* Change the text `ADDHOSTNAMEHERE` in the copy-pasted code to your machine's hostname.
* Copy and paste the address of the **source directory** after the text `export ATLASLC_SOURCEDIR=`.
* Copy and paste the address of the **data directory** after the text `export ATLASLC_DATA==`.
* Save the `atlaslc.sourceme` file.
* In the **data directory**, create a text file called `snlist.txt` that will house your SN list. Another option is to download the example file that is already set up with data and move it to the data directory.
* In the **source directory**, open `precursor.cfg`. This is the configuration file.
* Optional: You can change the `outsubdir` name in this file if you want to create different sub-directories for different SNe.
* After `username`, add your ATLAS machine username.
* Change any other values you wish (for example, the filter to be used).

### Usage (NEEDS UPDATING):
This code allows you to download SN light curves. There are many options to configure. You can view the list of arguments in the `SNloop.py` file in the source directory. To summarize:
* `download_lc_loop.py` initializes the program.
* Add the SN name(s) that you want light curves for (or use 'all' if you want to use all SNe in the `snlist.txt` file) to the command.
* `--passwd` specifies your password for the ATLAS machines (**required**). It is recommended to put your password in `'quotation marks'` to avoid any errors with symbols.
* `-s` will make sure all files you download are saved (**recommended; should be added to the command every time you run the code**).
* `-v` specifies verbose level (**recommended**).
* `-l` specifies the lookback time in days.
* `--api` lets you download light curves using the ATLAS API instead of SSH. You may have to change your username in the `precursor.cfg` file to match your API credentials.  
* `--snlistfilename` gives you the option to read a different text file that houses a SN list instead of using the address in the `precursor.cfg` file.
* `--user` specifies your username. If you use this argument, it will override the username you set in `<precursor.cfg>`.
* `--outrootdir` and `--outsubdir` give you the option to specify different output directories and sub-directories instead of the ones you specified in the `atlaslc.sourceme` file.
* `-c` lets you reference a different configuration file instead of the `precursor.cfg` file.
* `-e` adds an additional configuration file.
* `-f` specifies the filter (o or c) to be used when downloading the light curves. By default, the code will reference `precursor.cfg`, but adding this argument will override that.

#### Additional usage:
You can also use forced photometry to get offset light curves, plot each SN's light curve, and average the light curve data. To do any of these, follow the following instructions.
* **To get forced photometry offsets automatically**, add `--forcedphot_offset True` to the command. In the `precursor.cfg` file, specify the `radii` to be used (you can add multiple) and the `n` number of offsets per radius.
* **To plot each SN's light curve automatically**, add `--plot True` to the command. You can also specify the x and y limits using the following arguments: `--xlim_lower`, `--xlim_upper`, `--ylim_lower`, and `--ylim_upper`.
* **To average the light curve data automatically**, add `--averagelc True` to the command. In the `precursor.cfg` file, specify the `MJDbinsize` to be used, OR add `--MJDbinsize` to the command.
* **To detect precursor bumps in the data**, add `--detectbumps True` to the command. You can also simulate bumps using Gaussian models by using the argument `--sim_gaussian` and specifying the MJDs to insert a bump, the peak magnitude, and the length. For example, `--sim_gaussian 58000,58325 19,20,20.5,21 30` would add simulated bumps at the MJDs 58000 and 58325. These bumps would have peak magnitudes of 19, 20, 20.5, and 21, all shown in different plots. They would also have lengths of 30.

Additional functionality enables you to do these tasks using existing data that has already been downloaded.
* **To plot each SN's light curve using existing data**, initialize the program using `plot_lc.py`, then add to the command the SN name(s) you want plotted. Again, you can also specify the x and y limits using the following arguments: `--xlim_lower`, `--xlim_upper`, `--ylim_lower`, and `--ylim_upper`.
* **To clean up existing SN data and flag in the Mask column**, run `cleanup_lc.py`, then add to the command the SN name(s) you want plotted.
* **To average the light curves using existing data**, initialize the program using `averageLC.py`, then add to the command the SN name(s) you want plotted. You can also override the MJDbinsize you set in `precursor.cfg` by adding `--MJDbinsize` to the command.

#### Using autoadd.py with your snlist.txt
Given a TNS name, `autoadd.py` can automatically retrieve the RA and Dec coordinates from the TNS and add the TNS name, RA, and Dec to `snlist.txt`. To run the script, follow the following instructions:
* If you have only the TNS name of your SN, use the following command and add the TNS name at the end of the command: `autoadd.py`.
* If you have the SN name (could be TNS name or other designation) AND the RA and Dec coordinates, initialize the program (`autoadd.py`) and add the SN name at the end of the command. Then, use the arguments `--ra` and `--dec` to specify the RA and Dec, and run the program.

#### Example commands:
* `autoadd.py 2020lse` adds the TNS name 2020lse, its RA, and its Dec to `snlist.txt`. Similarly, `autoadd.py 2020lse --ra 10:41:02.20 --dec -27:05:00.3` adds the TNS name 2020lse, the given RA, and the given Dec to `snlist.txt`.
* `download_lc_loop.py 2020lse -v -s -l 70 --forcedphot_offset True --plot True --averagelc True --MJDbinsize 20 --passwd 'XXX'` gets the data for SN 2020lse with verbose level 1, overwrites files with the same name, saves the files, and uses a lookback time of 70. Then the offset data is downloaded, plots are saved, and the SN light curve data is averaged with an MJDbinsize of 20.
* `plot_lc.py 2020lse -v -s` plots the data for SN 2020lse as long as there is already existing light curve data.
* `offsetstats.py 2020lse -v -s --MJDbinsize 20` averages the data for SN 2020lse as long as there is already existing light curve data.
