# SAFETRANS visibility module (docker) readme

Description

The dockerized tool is set to make use of a folder containing a set of LIDAR measurements in order to calculate visibility ranges. File format must be as directed by RAYMETRICS manuals for LIDAR scan measuremets, with a 4-line metadata header.

## Usage:
Docker creation (requires "Docker"):
> git clone https://github.com/Hector-Xavier/SAFETRANS_LIDAR_analysis

> docker build -t safetrans SAFETRANS_LIDAR_analysis

Container deployment:
**Deployment example**
> docker run --rm -v /your/LIDAR/scan/folder:/data safetrans Rscript SAFETRANS_visibility_module.R 1 FALSE azimuth null FALSE null null TRUE TRUE

**"Silent" deployment example**
> docker run --rm -v /your/LIDAR/scan/folder:/data safetrans Rscript SAFETRANS_visibility_module.R 1 FALSE azimuth null FALSE null null TRUE TRUE > /dev/null

**Explicit form of arguments**
> docker run --rm -v /your/LIDAR/scan/folder:/data safetrans Rscript SAFETRANS_visibility_module.R channels is_scan scan_type model incoming height distance output verbose

**Granting user access to newly-created files**
> sudo chown -R $(id -u) your/LIDAR/scan/folder/Output


## Arguments:
**your/LIDAR/scan/folder** - Full path to the folder containing the LIDAR measurements.

**/data** - Do not chage: internal container directory interfacing with the external environment.

**safetrans** - The docker image ID, or its more human-readable _tag_, specified during creation.

**Rscript SAFETRANS_visibility_module.R** - Do not change: points to the script performing the analysis.

**channels** - Designates the channel or channels of the data scan files that will be used for analysis. At present, elevation scans should be made using only a single channel.

_Format examples:_
> 1

> 3

>3_4

> 1_2_3_4

**is_scan** - Designates whether the files present in the folder should be considered as parts of a scan profile, either azimuth or elevation. If set to elevation, the cartesian extinction profile is calculated and visibility ranges of incoming and outcoming objects at set positions for that azimuth angle become available.

_Format:_
> TRUE

> FALSE

**scan_type** - Designates the type of scan if _is_scan_ is set to TRUE.

_Format:_
> azimuth

> elevation

**model** - Designates the atmospheric model used for the translation of extinction coefficients from 1054 nm to the visible part of the spectrum. If set to null, no conversion is carried out and the results correspond to visibility ranges for 1054 nm.

_Format:_
> null

> maritime

> urban-rural

**incoming** - If _is_scan_ is set to TRUE and _model_ is set to elevation, it designates whether the object of interest is incoming (TRUE) or outcoming (FALSE).

_Format:_
> TRUE

> FALSE

**height** - If _is_scan_ is set to TRUE and _model_ is set to elevation, it designates the height of the object of interest (in metres).

_Format examle:_
> 4500

**distance** - If _is_scan_ is set to TRUE and _model_ is set to elevation, it designates the distance of the object of interest (in metres).

_Format example:_
> 3000

**output** - If set to TRUE, results are saved in output files (see next section).

_Format:_
> TRUE

> FALSE

**verbose** - If set to TRUE, various messages pertaining to the analysis process are displayed.

_Format_
> TRUE

> FALSE


## Output files

**Radial_extinction_coefficients_1054_nm.txt** - Tab-delimited text file containing the calculated extinction coefficients for all channels and all files designated as part of the analysis. If repeaded queries are made concerning the same initial data, it is read to reduce answer calculation delay.

**Radial outward visibility distance** - Teb-delimited text file containing the calculated radial visibility ranges for all channels and all files designated as part of the analysis. It is specific to the atmospheric _model_ used.

**Cartesian extinction profile** - Tab-delimited text file containing the extension of the radial extinction coefficient profile to cartesial coordinates. Due to its potentially great size, it is gz-compressed upon creation. **Only available for elevation scans.** It is specific to the atmospheric _model_ used.

**Incoming/Outcoming object visibility** - Tab-delimited text file containing visibility ranges of an incoming or outcoming object (respectively) at a set height and distance from the LIDAR location. Incoming/outcoming status, height and distance are supplied by the user and the file is overwritten after each query. **Only available for elevation scans.** It is specific to the atmospheric _model_ used.
