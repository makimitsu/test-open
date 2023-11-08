xCC = [50
60
100
150
200
300
400
500
600
700
800
850
1000
1200
1300
1500
2500
3500
4500
5500];
yCC = [0.0001
0.0005
0.0042
0.0197
0.0334
0.0749
0.1062
0.1284
0.1429
0.1526
0.1702
0.1806
0.196
0.2143
0.2175
0.2346
0.274
0.2761
0.2842
0.2817];
zCC=[0.0298
0.0295
0.023
0.0215
0.0197
0.0177
0.0161
0.0123
0.0127
0.0121
0.0109
0.0098
0.0093
0.0088
0.0087
0.0075
0.0063
0.0046
0.0046
0.0036];
xDC= [30
40
50
60
75
100
125
150
175
200
225
250
400
800
1000
1500
];
yDC=[0.0002
0.001
0.0014
0.0043
0.0057
0.0102
0.0119
0.013
0.0176
0.0161
0.0215
0.0208
0.0195
0.0196
0.0171
0.0161
];
zDC=[0.3159
0.2998
0.2873
0.2781
0.2621
0.2413
0.2393
0.2278
0.2191
0.2146
0.1941
0.1964
0.1755
0.1313
0.1175
0.0962
];
xCW = [2000
1500
1000
900
800
700
600
500
400
300
200
175
150
125
100
75
60
];
yCW=[0.44789
0.38974
0.30856
0.289
0.26676
0.24307
0.21275
0.18201
0.14628
0.10397
0.05486
0.04334
0.03311
0.01962
0.0099
0.00212
0.0001];
zCW=[0.441
0.453
0.476
0.482
0.487
0.494
0.505
0.513
0.528
0.539
0.564
0.571
0.58
0.592
0.606
0.628
0.644];
xDW=[300
400
500
600
700
800
900
1000
1200
1500];
yDW=[0.0001
0.0013
0.0014
0.0028
0.0043
0.0037
0.0055
0.0052
0.0072
0.0082];
zDW=[0.6162
0.6005
0.5921
0.5823
0.5638
0.5584
0.5576
0.5479
0.5237
0.508];

xNeC=[70
100
125
150
200
300
400
500
600
700
850
1000
1200];
yNeC=[0.0001
0.0023
0.0056
0.0118
0.0271
0.0631
0.1028
0.1315
0.1563
0.1844
0.2115
0.2468
0.281];
zNeC=[0.0021
0.0016
0.0012
0.0016
0.0014
0.0014
0.0019
0.0014
0.0016
0.0014
0.0013
0.0012
0.0008];
xNeW=[600
500
400
300
200
150
100
75
50
40
];
yNeW=[0.4325
0.3678
0.293
0.21
0.1206
0.07801
0.03098
0.01159
0.00098
0.00006];
zNeW=[0.463
0.472
0.488
0.504
0.529
0.549
0.582
0.603
0.644
0.662];
xHeC=[30
40
50
60
100
200
400
800];
yHeC=[0.0001
0.0015
0.0033
0.0052
0.0179
0.0412
0.0585
0.0677];
zHeC=[0.2116
0.197
0.1946
0.1889
0.1689
0.1291
0.1071
0.0816];
xHeW = [
144
145
150
175
200
300
400
500
600
700
800
900
1000];
yHeW=[0.0001
0.0002
0.0003
0.0012
0.0013
0.0075
0.0162
0.0215
0.027
0.0318
0.0372
0.0387
0.0434];
zHeW=[0.6648
0.6598
0.6549
0.6456
0.6361
0.6215
0.611
0.6007
0.5838
0.5825
0.5755
0.5722
0.5768];

lw = 1.5;
figure
% loglog(xCC,yCC,'r+--','LineWidth',lw)
% hold on
loglog(xDC,yDC,'b+--','LineWidth',lw)
hold on
loglog(xNeC,yNeC,'g+--','LineWidth',lw)
hold on
loglog(xHeC,yHeC,'k+--','LineWidth',lw)
hold on
loglog(xCW,yCW,'ro-','LineWidth',lw)
hold on
loglog(xDW,yDW,'bo-','LineWidth',lw)
hold on
loglog(xNeW,yNeW,'go-','LineWidth',lw)
hold on
loglog(xHeW,yHeW,'ko-','LineWidth',lw)
xlim([10 10000])
ylim([1E-4 1])
legend('D on C','Ne on C','He on C','C on W','D on W','Ne on W','He on W')
xlabel('Incident Energy [eV]')
ylabel('Sputtered Yield')
fontsize(16,"points")
grid on
hold off

figure
semilogx(xCC,zCC,'r+--','LineWidth',lw)
hold on
semilogx(xDC,zDC,'b+--','LineWidth',lw)
hold on
semilogx(xNeC,zNeC,'g+--','LineWidth',lw)
hold on
semilogx(xHeC,zHeC,'k+--','LineWidth',lw)
hold on
semilogx(xCW,zCW,'ro-','LineWidth',lw)
hold on
semilogx(xDW,zDW,'bo-','LineWidth',lw)
hold on
semilogx(xNeW,zNeW,'go-','LineWidth',lw)
hold on
semilogx(xHeW,zHeW,'ko-','LineWidth',lw)
legend('C on C','D on C','Ne on C','He on C','C on W','D on W','Ne on W','He on W')
xlabel('Incident Energy [eV]')
ylabel('Reflection Ratio')
fontsize(16,"points")