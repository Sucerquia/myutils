%chk=G-stretched08
%NProcShared=8
#P bmk/6-31+g ! ASE formatted method and basis
opt(modredun,calcfc)

Gaussian input prepared by ASE

0 1
C                 0.0000000000        0.0000000000        0.0000000000
H                -0.1266115579        1.0855567452        0.0001249573
H                -0.3920554370       -0.4720084529       -0.9060521324
H                -0.3920750844       -0.4722214196        0.9059318676
C                 1.7514794087       -0.5285542067       -0.0000374087
O                 1.8168008568       -1.7663837538       -0.0001253866
N                 3.0321520303        0.1971686305        0.0000120237
H                 3.0422809249        1.2095799140        0.0000890271
C                 4.4688206831       -0.5168425666       -0.0000554912
H                 4.4857522677       -1.1639102448       -0.8862664855
H                 4.4857427650       -1.1640920516        0.8860215145
C                 5.9549748162        0.3464564822        0.0000380105
O                 5.9089975905        1.5786334443        0.0001709950
N                 7.2245016029       -0.3825515214       -0.0000415609
H                 7.1111428478       -1.3907356292       -0.0001545991
C                 8.8856973998        0.0000000000        0.0000000000
H                 9.2986528966       -0.4392410387       -0.9131628606
H                 8.9490914219        1.0895592025        0.0001000214
H                 9.2986480346       -0.4394083109        0.9130851394

1 16 F
