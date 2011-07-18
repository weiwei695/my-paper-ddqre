
%--------------变调不变语速程序---------------------------%
function [xfinal,Fs]=Pitch_Modification(x,f,semitone)
% xfinal:最终TSM后的数据
% x:需要TSM的数据
% scale：时间缩放比率,大于1变快，小于1变慢 

Fs=f;

semitone = round(semitone);
if semitone>12||semitone<-12
    semitone = 0;
    warndlg('semitone has to be in [-12 12],semitone is reset to be 0');
end

sent_L = length(x);
scale = 2^(semitone/12);
L = 160; %frame length
S = L/4; %hop size

overlap = L - S;
Nframe = floor((sent_L-overlap)/S);
Lq = round(L*scale);

win=hamming(L);%或者用下面5条语句，也即构造hamming window
% a = 0.54;
% b = -0.46;
% n = 1:L;
% win = sqrt(S)/sqrt((4*a^2+2*b^2)*L)*(a+b*cos(2*pi*n/L));
% win = win(:);

winq=hamming(Lq);%或者用下面3条语句，也即构造hamming window
% n = 1:Lq;
% winq = sqrt(S)/sqrt((4*a^2+2*b^2)*Lq)*(a+b*cos(2*pi*n/Lq));
% winq = winq(:);

Nit = 4;%迭代次数using in 函数iterated_recon中
xfinal  = zeros(sent_L,1);
U = sum(win)/(S);

k = 1;
% kk = 1;
h = waitbar(0,'Please wait...');
for n = 1:Nframe
    if k:k+Lq-1<=sent_L
        frm = winq.*x(k:k+Lq-1)/U;
    else
        frm = winq.*[x(k:sent_L);zeros(Lq - (sent_L-k+1),1)]/U;
    end
    
    if 1
        a = lpc(frm,32);
        est_frm = filter([0 -a(2:end)],1,frm);
        e = frm - est_frm;
        frm = e;
    end

    frm_resamp = resample(frm, L, Lq);
    xSTFTM = abs(fft(frm_resamp));

    if k+L-1<=sent_L
        res = xfinal(k:k+L-1);
    else
        res = [xfinal(k:sent_L);zeros(L - (sent_L-k+1),1)];
    end

    x_recon = iterated_recon(xSTFTM, res, Nit, win);
    
    
    if 1
        x_recon = filter(1,a,x_recon);
    end

    
    if (k+L-1<=sent_L)
        xfinal(k:k+L-1) = xfinal(k:k+L-1) + x_recon;
    else
        xfinal(k:sent_L) = xfinal(k:sent_L) + x_recon(1:sent_L-k+1);
    end
    k = k + S;

    waitbar(n/Nframe, h)
end
close(h);

function x_recon = iterated_recon(xSTFTM, x_res, Nit, win)
j = sqrt(-1);
for i = 1:Nit
    phi = phase(fft(win.*x_res)) + randn(size(x_res))*0.01*pi;% random phase purturbation will reduce some resonance.
    x = xSTFTM.*exp(j*phi); %M-constraint
    x_recon = ifft(x);
    x_res = real(x_recon);
end
x_recon = x_res;


