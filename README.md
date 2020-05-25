# ExoVis

This repository contains the code of a prototype for our paper on progresive parameter space visualization.
It is published as supplementary material to our paper _"Progressive Parameter Space Visualization for Task-Driven SAX Configuration"_, which will be presented at this year's [EuroVA workshop](https://www.eurova.org/eurova-2020).

The focus of this work is a parameter space matrix, which progressively generates an error estimation of computing a time series analytic on data aggregated using [SAX](https://www.cs.ucr.edu/~eamonn/SAX.htm) and [PAA](https://jmotif.github.io/sax-vsm_site/morea/algorithm/PAA.html).
Each cell of this matrix represents the error of using a particular combination of a SAX alphabet size ```alpha``` and a PAA division integer ```omega```, compared to running the same computation on a sample of the unaggregated data.
This allows you to make an informed decision, as to which parameter combination produces useful results while also significantly reducing computation times.
This error is gradually updated, as more and more data is processed for that computation.
This progression is indicated by progress bars along the columns and rows of the matrix and through a tooltip that shows the per-cell value.

You may also steer the computation towards regions of interest by either clicking an individual cell or by selecting a full column or row by clicking on the respective progression bar.
Then, these cells are fully processed first, which can reduce the waiting times of computing the full matrix.


## Getting started
In order to run the prototype, a Python environment with version 3.7 or higher is needed.
We recommend using a basic [Anaconda installation](https://docs.conda.io/en/latest/miniconda.html), which is available for all platforms.

### Installing dependencies

Then, we need to install the required dependencies.
We do so by running in the terminal:

```bash
pip install lightkurve numpy astropy eel saxpy scipy pyDOE
```


### Launching the backend

This repository comes with two example use cases for calculating the computation error on the aggregated data.
The first uses the light intensity measured by NASA during the [Kepler mission](https://www.nasa.gov/mission_pages/kepler/main/index.html) in the search for exoplanets.
Here, the light curve data is scanned for periodic dips, that indicate a planet transiting a star in regular intervals.

To launch the script for this use case run the following command in the ```backend/``` directory of this repository:
```
python server_lightcurve.py
```

The second use case analyses weather data from different countries for their general climate by fitting a regression line to the temperature data measured over a period of around 200 years.

To launch the script for this use case run the following command in the ```backend/``` directory of this repository:
```
python server_weather.py
```

When running these scripts for the first time may take a while, as the required data will be downloaded first.
Once the data is downloaded, the terminal will show the following message:
```
Backend launched successfully. Waiting for requests ...
```

### Launching the frontend

Now, launch the frontend in your browser.
For this, you will need to open one of the two .html files in the ```frontend/``` directory of this repository.
Depending on which backend script you have launched in your terminal, open either ```index_lightcurve.html``` or ```index_weather.html``` in your browser of choice.
The matrix should now start generating, with the backend script in the terminal notifying you about its progress.


## Citing this work
In order to cite this project in your work, we recommend using the following bibtex entry:

```bib
@inproceedings {Loeschcke,
  booktitle = {EuroVis Workshop on Visual Analytics (EuroVA)},
  editor = {Turkay, Cagatay and Vrotsou, Katerina},
  title = {{Progressive Parameter Space Visualization for Task-Driven SAX Configuration}},
  author = {Loeschcke, Sebastian and Hogr\"{a}fer, Marius and Schulz, Hans-J\"{o}rg},
  year = {2020},
  publisher = {The Eurographics Association},
  ISSN = {2664-4487},
  ISBN = {978-3-03868-116-8},
  DOI = {10.2312/eurova.20201085}
}
```



