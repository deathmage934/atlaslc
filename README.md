# PROJECT NAME: atlaslc

## Description: 
Using a table of SN names, RA, and Dec, atlaslc logs into the ATLAS machines and retrieves the SN light curve data. You will need a table containing your **SN names, RA and Dec for each SN, and your username and password** for the ATLAS machines.

## Installation:
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

## Usage:
This code allows you to download SN light curves. There are many options to configure. You can view the list of arguments in the `SNloop.py` file in the source directory. To summarize:
* `download_lc_loop.py` initializes the program.
* Add the SN name(s) that you want light curves for (or use 'all' if you want to use all SNe in the `snlist.txt` file) to the command.
* `--snlistfilename` gives you the option to read a different text file that houses a SN list instead of using the address in the `precursor.cfg` file.
* `-s` will make sure all files you download are saved (**recommended; should be added to the command every time you run the code**).
* `-v` specifies verbose level (**recommended**).
* `-l` specifies the lookback time in days.
* `-o` will overwrite any files with the same name (**recommended** if you download data for the same SN multiple times).
* `--user` specifies your username. If you use this argument, it will override the username you set in `<precursor.cfg>`.
* `--passwd` specifies your password for the ATLAS machines (**required**). It is recommended to put your password in 'quotation marks' to avoid any errors with symbols.
* `--outrootdir` and `--outsubdir` give you the option to specify different output directories and sub-directories instead of the ones you specified in the `atlaslc.sourceme` file.
* `-c` lets you reference a different configuration file instead of the `precursor.cfg` file.
* `-e` adds an additional configuration file.
* `-f` specifies the filter (o or c) to be used when downloading the light curves. By default, the code will reference `precursor.cfg`, but adding this argument will override that.

### Additional usage:
You can also use forced photometry to get offset light curves, plot each SN's light curve, and average the light curve data. To do any of these, follow the following instructions.
To get forced photometry offsets automatically:
* Add `--forcedphot_offset True` to the command.
* In the `precursor.cfg` file, specify the `radii` to be used (you can add multiple) and the `n` number of offsets per radius.
To plot each SN's light curve automatically:
* Add `--plot True` to the command.
To average the light curve data:
* Add `--averagelc True` to the command.
* In the `precursor.cfg` file, specify the `MJDbinsize` to be used.
Additional functionality enables you to do these tasks using existing data that has already been downloaded.
To plot each SN's light curve using existing data:
* `plot_lc.py` initializes the program. Add to the command the SN name(s) you want plotted.
To average the light curves using existing data:
* `averagelc_loop.py` initializes the program. Add to the command the SN name(s) you want plotted.
* Optional: you can also override the MJDbinsize you set in `precursor.cfg` by adding `--MJDbinsize` to the command.

### Example commands:
* `download_lc_loop.py 2020lse -v -o -s -l 70 --forcedphot_offset True --plot True --averagelc True --MJDbinsize 20 --password 'XXX'` gets the data for SN 2020lse with verbose level 1, overwrites files with the same name, saves the files, and uses a lookback time of 70. Then the offset data is downloaded, plots are saved, and the SN light curve data is averaged with an MJDbinsize of 20.
* `plot_lc.py 2020lse -v -s` plots the data for SN 2020lse as long as there is already existing light curve data.
* `averagelc_loop.py 2020lse -v --MJDbinsize 20` averages the data for SN 2020lse as long as there is already existing light curve data.
