## **PSCF-CWT Dependencies:**
  •	Pandas  
  •	Numpy  
  •	Sys  
  •	Os  
  •	Math  
  •	Datetime
## Data Preparation:  
•	Prepare the HYSPLIT tdump file.  
•	Prepare pollutant concentration file, the time data should be in the first column, and the pollutant concentration data should be in the second column.
## Input File Paths:  
•	On lines 321-322, enter the paths for the tdump folder and the pollutant file.  
•	On lines 323-324, enter the output file paths for the HYSPLIT trajectory to pollutant matrix conversion and the CWT calculation.
## Configuration:  
•	Modify the grid point resolution and weight function values on lines 331 and 332, respectively:  
•	3n_ave < n_ij  
•	1.5n_ave < n_ij ≤ 3n_ave  
•	n_ave < n_ij ≤ 1.5n_ave  
•	n_ij < n_ave
## Running the Code:  
•	Call the File_load_HYSPLIT function to align the tdump file with the pollutant file based on time, and convert the data into a matrix corresponding to latitude, longitude, and grid resolution. Store the matrix in the pwd_matrix_output folder.  
•	Call the CWT function to perform the CWT calculation and store the results in the pwd_matrix_CWT_output folder.
