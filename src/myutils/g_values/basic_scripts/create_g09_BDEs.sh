python ../create_g09_file.py opt.xyz g09BDE

sed -i "1a %NProcShared=8" g09BDE.com
sed -i "s/^0 1$/0 2/g" g09BDE.com
sed -i "s/%chk=test/%chk=g09BDE/g" g09BDE.com
sed -i "3a opt(modredun,calcfc) freq" g09BDE.com
