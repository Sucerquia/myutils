%chk=G-stretched09
%NProcShared=8
#P bmk/6-31+g ! ASE formatted method and basis
opt(modredun,calcfc)

Gaussian input prepared by ASE

0 1
C                 0.0000000000        0.0000000000        0.0000000000
H                -0.1899088695       -1.0704434897       -0.0000238746
H                -0.1269396191        0.5492523854        0.9284810838
H                -0.1269392377        0.5492924031       -0.9284589162
C                 2.4093077619        0.5223730536        0.0000644089
O                 2.3653807236        1.7464912939        0.0002174380
N                 3.6471796709       -0.1608145362       -0.0000214085
H                 3.6752766738       -1.1719341178       -0.0001494271
C                 4.9930088700        0.5817012932        0.0000617027
H                 5.0066336585        1.2319584408        0.8859476937
H                 5.0066369088        1.2321518276       -0.8856823063
C                 6.3648887356       -0.2826864658       -0.0000352032
O                 6.3182847952       -1.5184472967       -0.0001821725
N                 7.5736867268        0.4455133599        0.0000549985
H                 7.4698280275        1.4532066363        0.0001790671
C                 9.0856975211       -0.0000000000        0.0000000000
H                 9.5508383417        0.3961123206        0.9063116928
H                 9.0998113848       -1.0930356167       -0.0001020093
H                 9.5508265480        0.3962839705       -0.9062423072

1 16 F