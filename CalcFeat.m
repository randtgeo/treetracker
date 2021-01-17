cloudpath = 'Z:\PersonalShare\Luke\BGC\';
file = 'MOBERLY_2020-05.txt'
outfile = 'MOBERLY_2020-05-new'
CCpath = 'C:\Program Files\CloudCompare';


cmd = ['CloudCompare -AUTO_SAVE OFF -O -GLOBAL_SHIFT AUTO ',cloudpath,file,...
    ' -REMOVE_ALL_SFS -REMOVE_RGB -REMOVE_NORMALS -FEATURE LINEARITY 1.0 -FEATURE LINEARITY 0.5 -FEATURE VERTICALITY 0.5' ...
    ' -DENSITY 0.25 -TYPE VOLUME -C_EXPORT_FMT ASC -SAVE_CLOUDS FILE ',cloudpath,outfile]

cd(CCpath)
system(cmd)
cd(cloudpath)