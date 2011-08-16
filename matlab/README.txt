= ABOUT MTKM =

This is the initial public release of the Manifold Toolkit for Matlab
(MTKM).  It accompanies the IROS 2011 paper:

  R. Wagner, O. Birbach, U. Frese: Rapid Development of Manifold-Based
  Graph Optimization Systems for Multi-Sensor Calibration and SLAM

Unless otherwise noted all code is made available under a 3-clause
BSD-style open source license.

If you have any trouble running the provided code or any other queries
please do not hesitate to contact Rene Wagner <rene.wagner@dfki.de>.


= RUNNING THE CALIBRATION EXAMPLES = 

The calibration examples are located under
  * examples/calibration_mono_codesample for the introductory
    calibration example (Section IV)
  * examples/calibration_justin_codesample for the multi-sensor cross
    calibration example (Section V)
  * examples/calibration_kinect for the Microsoft Kinect calibration
  * examples/calibratiob_nao for the NAO calibration

For executing each of the individual examples, start Matlab, cd to the
desired folder, run the included setuppath script in that
folder and run the script by calling one of

  * monocalib
  * justincalib
  * kinectcalib
  * naocalib

respectively. 


= RUNNING THE SLAM EXAMPLES =

The SLAM examples are located under

  * examples/spa for the SPA data set example
  * examples/dlr-spatial_cognition for the DLR Spatial Cognition data
    set example
  * examples/toro3d for the synthetic sphere data set example
  * examples/g2o3D for the parking garage data set example

Before running any of the SLAM examples please download the respective
data set as indicated by the DATASET.txt file in each example folder.

Then start Matlab, cd to one of the example folders, run the
setuppaths script in that folder and run the example script, i.e. one
of

  * spa_test
  * dlr_test
  * toro3d_test
  * g2o_test

Each script will parse the respective data set, create the
optimization problem, run the optimizer (depending on the example
either Levenberg-Marquardt or Gauss-Newton), and plot the resulting
map after optimization.

If the parsing fails, please make sure you have downloaded the
respective data set to the current directory and that the file name
matches the one expected by the *_test.m script.

Please note that the SLAM examples require a lot of memory.  We found
that 4GB of main memory (RAM) is sufficient for all examples.  Make
sure that the kernel was compiled with PAE support when working on a
32 bit Linux machine.  Many Linux distributions provide such kernels
as part of separate -PAE kernel packages.


= REPRODUCING THE SLAM PLOTS =

You will need the exportfig toolbox to generate the EPS files used in
the paper

  http://www.mathworks.com/matlabcentral/fileexchange/727

Download and extract it and add the directory to your Matlab
path. Then, run the following scripts:

  * Fig. 1 (parking garage data set):
  
      examples/g2o3D/g2o_plots.m

  * Fig. 6 (SPA samples):
     
      examples/spa/spa_plots.m

  * Fig. 7 (DLR Spatial Cognition data set):

      examples/dlr-spatial_cognition/dlr_plot_before_after.m

  * Fig. 8 (synthetic sphere):

      examples/toro3d/plot_sphere_before_after.m

All EPS files are dumped in the current directory. Also see the note
on the setuppaths scripts above.
