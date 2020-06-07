An automated method of quantifying particle velocity from linescans is very
useful for analyzing capillary-level blood flow data obtained with two-photon
microscopy.  The Radon transform method is an fast and robust method for
obtaining velocity information from linescan data. 

This demo code generates pseudo-data, and then uses the Radon transform to
calculate the angle of the streaks in a user-defined time window.  When the
spatial dimension of the space-time plot is along the x axis, and the time
dimension is along the y-axis,  the cotangent of the angle of the streaks is
proportional to the velocity.  The error as a fraction of the velocity will be
largest for near horizontal and vertical lines.  For angles of
less than 10 degrees from horizontal, a cross-correlation method is more
appropriate.  


The publication describing the method is described here:

-   Drew PJ, Blinder P, Cauwenberghs G, Shih AY, Kleinfeld D, Rapid
    determination of particle velocity from space-time line-scan data using the
    Radon transform, Journal of Computational Neuroscience, 29(1-2):5-11

    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4962871/>

 

 

See also:

-   Chhatbar PY, Kara P. Improved blood velocity measurements with a hybrid
    image filtering and iterative Radon transform algorithm. *Frontiers in
    Neuroscience*. 2013;7:106. doi:10.3389/fnins.2013.00106.
