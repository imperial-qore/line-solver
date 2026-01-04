function model = jsimgEdit
fname = [lineTempName,'.jsimg'];
copyfile(which('emptyModel.jsimg'),fname);
jsimgView(fname);
try
    model = JMT2LINE(fname);
    delete(fname);
catch
    line_error(mfilename,'Error: JMT2LINE failed. JMT model stored at: %s\n',fname)
end
end