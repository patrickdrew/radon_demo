
function []=radon_demo(~)
% radon demo This code generates pseudo-data to mimic the data obtained
% with line-scans along vessels and then uses the Radon transform to find
% the angle of the 'streaks'.  The pseudo-data mimics a wide variety of
% flow 'velocities'.

% The Radon method is described in detail in:
% Drew PJ, Blinder P, Cauwenberghs G, Shih AY, Kleinfeld D, Rapid
% determination of particle velocity from space-time line-scan data using
% the Radon transform, Journal of Computational Neuroscience, 29(1-2):5-11
% https://sites.esm.psu.edu/~pjd17/Drew_Lab/Resources.html
% Patrick Drew, pjd17@psu.edu
% 

% code for generating simulated data, parameters chosen to generate wide variety of 'velocities' 
spatial_freq=1;
npoints=64*spatial_freq;%number of points in a line
nlines=64*8*2*spatial_freq;%number of lines 
zz=zeros(nlines,npoints);
f0=30;
phi0=(pi/25);
phi1=1;
f_phi=pi;
for i=1:npoints
    for j=1:nlines
        phi=phi0*i;%
        f=f0+phi1*sin(((f_phi*2*pi*j)/spatial_freq)/1000);
        zz(j,i)=0.5+0.5*sin((2*pi*f*j/spatial_freq)/1000+phi);
    end
end

%plots the artificial 'linescan'
close all
figure(1)
ax(1)=subplot(211);
imagesc(zz')
colormap gray
xlabel('time, A.U.')
ylabel('space, A.U.')

[thetasz32,the_tz32,spread_radon32]=GetVelocityRadonFig_demo(zz,npoints);
ax(2)=subplot(212);
plot(the_tz32,thetasz32,'b.-')% plots the angle at any given time point
linkaxes(ax,'x')
hold on
xlabel('time, A.U.')
ylabel('angle, degrees')


function [thetas,the_t,spread_matrix]=GetVelocityRadonFig_demo(data,windowsize);
%calculates the angle of the streaks in a space-time 'data' within a time window set by windowsize
%after an initial rough calulation of the angle, the angle is recalulated
%around the rough peka with
%OUTPUTS
%thetas - the time varying agle of the space-time image
%the_t - time pointsof the angle estimates (in lines)
%spreadmatrix - matix of variances as a function of angles for each window
%position
%INPUTS
%data - the matrix of time X space data 
%windowsize - number of lines to use in estimating velocity.  must be a
%multiple of 4
stepsize=.25*windowsize;
nlines=size(data,1);
npoints=size(data,2);
nsteps=floor(nlines/stepsize)-3;
%find the edges
angles=(0:179);
angles_fine=-2:.25:2;

spread_matrix=zeros(nsteps,length(angles));
spread_matrix_fine=zeros(nsteps,length(angles_fine));
thetas=zeros(nsteps,1);

hold_matrix=ones(windowsize,npoints);
the_t=NaN*ones(nsteps,1);

for k=1:nsteps
    the_t(k)=1+(k-1)*stepsize+windowsize/2;
    data_hold=data(1+(k-1)*stepsize:(k-1)*stepsize+windowsize,:);
    data_hold=data_hold-mean(data_hold(:))*hold_matrix;%subtract the mean
    radon_hold=radon(data_hold,angles);%radon transform
    spread_matrix(k,:)=var(radon_hold);%take variance
    [~, the_theta]=max(spread_matrix(k,:));%find max variace
    thetas(k)=angles(the_theta);     
    radon_hold_fine=radon(data_hold,thetas(k)+angles_fine);%re-do radon with finer increments around first estimate of the maximum
    spread_matrix_fine(k,:)=var(radon_hold_fine);
    [~, the_theta]=max(spread_matrix_fine(k,:));
    thetas(k)=thetas(k)+angles_fine(the_theta);
end
thetas=thetas-90; %rotate

