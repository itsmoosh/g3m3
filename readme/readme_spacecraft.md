Spacecraft may optionally be included in multifluid simulations. The presence of spacecraft is toggled on in the input file via the 'spacecraft' boolean flag.

Spacecraft position/time values are read in from files in the multifluid/spacecraft_info directory. Data recorded by spacecraft are stored in multifluid/data/spacecraft_data/. An output recording is generated for each position/time file in spacecraft_info/.

In addition, several default spacecraft are positioned at fixed locations relative to the parent body. These locations are specified in multifluid/spacecraft_info/defaults.pos .
