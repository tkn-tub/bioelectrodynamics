FIVE-REACTION DIFFERENTIATORS: OCTOBER 6 KINETICS DESIGN
=====================================================

Install and rerun
-----------------

For looking in Brenda database https://www.brenda-enzymes.info/
Place the file on your MATLAB path. With a supported Python installation configured through pyenv, install Zeep once:
pe = pyenv;
assert(strlength(pe.Executable)>0, 'Configure Python using pyenv first.');

[status,msg] = system(sprintf('"%s" -m pip install zeep', ...
                             char(pe.Executable)));
assert(status==0, msg);

You also need to input your credentials for log-in in Brenda website in the file BRENDA_credentials.txt with your email and password.

