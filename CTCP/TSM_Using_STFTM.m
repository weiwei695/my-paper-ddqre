%----变速不变频程序算法,采用Short Time Fourier Transform Magnitude Spectra---%
function [xfinal,Fs]=TSM_Using_STFTM(x,f,scale)
% xfinal:最终TSM后的数据
% x:需要TSM的数据
% scale：时间缩放比率,大于1变快，小于1变慢 

Fs=f;

if(scale>2||scale<0.5)
    scale = 1;
    warndlg('scale must be limited in [0.5 2],scale is reset to be 1');
end

sent_L=length(x);
L=256; %frame length
S=L/4; %hop size
m_S=round(S/scale);
overlap=L-S;
Nframe=floor((sent_L-overlap)/S);
Nit=5; %迭代次数using in 函数iterated_recon中

win=hamming(L);%或者用下面5条语句，也即构造hamming window
% a = 0.54;
% b = -0.46;
% n  = 1:L;
% win = sqrt(S)/sqrt((4*a^2+2*b^2)*L)*(a+b*cos(2*pi*n/L));
% win  = win(:);

L_recon=round(sent_L/scale);%变速后语音数据长度
xfinal=zeros(L_recon,1);% xfinal存贮变速后语音数据
U=sum(win)/(m_S);


k=1;
kk=1;
h=waitbar(0,'Please wait...');
for n = 1:Nframe
    frm = win.*x(k:k+L-1)/U;
    
    
    if 1
        a = lpc(frm,12);
        est_frm = filter([0 -a(2:end)],1,frm);
        e = frm - est_frm;
        frm = e;
    end
    
    
    xSTFTM = abs(fft(frm));
    if(kk+L-1<=L_recon)
        res = xfinal(kk:kk+L-1);
    else
        res = [xfinal(kk:L_recon);zeros(L - (L_recon-kk+1),1)];
    end
    x_recon = iterated_recon(xSTFTM, res, Nit, win);
    
    
    if 1
        x_recon = filter(1,a,x_recon);
    end
    
    
    if (kk+L-1<=L_recon)
        xfinal(kk:kk+L-1) = xfinal(kk:kk+L-1) + x_recon;
    else
        xfinal(kk:L_recon) = xfinal(kk:L_recon) + x_recon(1:L_recon-kk+1);
    end
    k = k + S;
    kk = kk + m_S;
    waitbar(n/Nframe, h)
end
close(h);


function x_recon = iterated_recon(xSTFTM, x_res, Nit, win)
j=sqrt(-1);
for i=1:Nit
    phi=phase(fft(win.*x_res)) + randn(size(x_res))*0.01*pi; % random phase purturbation will reduce some resonance.
    x=xSTFTM.*exp(j*phi);
    x_recon=ifft(x);
    x_res=real(x_recon);
end
x_recon = x_res;






