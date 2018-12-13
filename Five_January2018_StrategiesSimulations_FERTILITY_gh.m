% makes Figure 3 panel 1

clear all
%(1)initialise here
nGen = 5000;
f4b = 4;
Ni = 100;

pv = 0.6; %0:0.01:0.9;%0.75;
ph = 0.6; %0:0.01:1;

s1 = 0.6;
s2 = 0.7;
s3 = 0.7;
s4 = 0.7;
s5 = 0.4;

w4i= 0:0.01:2;% fertility benefit
ws = 0;% survival benefit

FREQ_STRAT1 = zeros(numel(w4i),5);
FREQ_STRAT2 = zeros(numel(w4i),5);

%FREQ_STRAT3 = zeros(numel(pvi),numel(phi),5);% hewlett 1
%FREQ_STRAT4 = zeros(numel(pvi),numel(phi),5);% hewlett 2
for m = 1:numel(w4i)
    w4 = w4i(m);
    %% Late horizontal
    %age class frequencies
    n1 = zeros(1,nGen);
    n2 = zeros(1,nGen);
    n3 = zeros(1,nGen);
    n4 = zeros(1,nGen);
    n5 = zeros(1,nGen);
    %informed individuals
    x1 = zeros(1,nGen);
    x2 = zeros(1,nGen);
    x3 = zeros(1,nGen);
    x4 = zeros(1,nGen);
    x5 = zeros(1,nGen);
    %first gen
    n1(1) = 0.2*Ni;
    n2(1) = 0.2*Ni;
    n3(1) = 0.2*Ni;
    n4(1) = 0.2*Ni;
    n5(1) = 0.2*Ni;
    x1(1) = 0.001;
    x2(1) = 0.001;
    x3(1) = 0.001;
    x4(1) = 0.001;% trait innovated by reproductive adults
    x5(1) = 0.001;
    
    n1(2) = 0.2*Ni;
    n2(2) = 0.2*Ni;
    n3(2) = 0.2*Ni;
    n4(2) = 0.2*Ni;
    n5(2) = 0.2*Ni;
    x1(2) = 0.001;
    x2(2) = 0.001;
    x3(2) = 0.001;
    x4(2) = 0.001;% trait innovated by reproductive adults
    x5(2) = 0.001;
    
    n1(3) = 0.2*Ni;
    n2(3) = 0.2*Ni;
    n3(3) = 0.2*Ni;
    n4(3) = 0.2*Ni;
    n5(3) = 0.2*Ni;
    x1(3) = 0.001;
    x2(3) = 0.001;
    x3(3) = 0.001;
    x4(3) = 0.001;% trait innovated by reproductive adults
    x5(3) = 0.001;
    
    finalPopSizeH_ph = zeros(5,1);
    finalXPropH_ph = zeros(5,1);
    finalAgeStructure_ph = zeros(5,1);
    
    Nh1000ph = zeros(1,1);
    
    N = Ni;
    Nh=zeros(1,1);
    
    %main loop vertical model
    
    for j = 4:nGen
        n1(j) = n4(j-1)*((f4b+w4)*x4(j-1)+f4b*(1-x4(j-1)));
        x1(j) = pv*((f4b+w4)*x4(j-1))/((f4b+w4)*x4(j-1)+f4b*(1-x4(j-1)));%these guys move to age class 2 now
        
        vert2 = (n4(j-2)*(s4+ws)^1*(f4b+w4)*x4(j-2))/(n4(j-2)*(s4+ws)^1*(f4b+w4)*x4(j-2)+n4(j-2)*f4b*s4^1*(1-x4(j-2)));% messed this up - FIX ME FIRST!!!
        n2(j) = (x1(j-1)*(s1+ws)+(1-x1(j-1))*s1)*n1(j-1);
        x2(j) = x1(j-1)+(1-x1(j-1))*(vert2*pv);
        
        vert3 = (n4(j-3)*(s4+ws)^2*(f4b+w4)*x4(j-3))/(n4(j-3)*(s4+ws)^2*(f4b+w4)*x4(j-3)+n4(j-3)*f4b*s4^2*(1-x4(j-3)));
        n3(j) = (x2(j-1)*(s2+ws)+(1-x2(j-1))*s2)*n2(j-1);
        x3(j) = x2(j-1)+(1-x2(j-1))*(vert3*pv);
        
        n4(j) = (x3(j-1)*(s3+ws)+(1-x3(j-1))*s3)*n3(j-1);
        x4(j) = x3(j-1)+(1-x3(j-1))*ph*((x3(j-1)*n3(j-1)+x4(j-1)*n4(j-1)+x5(j-1)*n5(j-1))/(n3(j-1)+n4(j-1)+n5(j-1)));
        
        n5(j) = n4(j-1)*x4(j-1)*(s4+ws)+n4(j-1)*(1-x4(j-1))*(s4)+n5(j-1)*x5(j-1)*(s5+ws)+n5(j-1)*(1-x5(j-1))*(s5);
        mixProp = (n4(j-1)*x4(j-1)*(s4+ws)+n4(j-1)*(1-x4(j-1))*(s4))/n5(j);
        x5(j) = x4(j-1)*mixProp+(1-mixProp)*x5(j-1);
        
        N = n1(j)+n2(j)+n3(j)+n4(j)+n5(j);
        
        if N<1
            break
            disp(['pv= ',num2str(pv),'and ph= '],num2str(ph))
        end
    end
    
    %x_save = [x1; x2; x3; x4; x5];
    %save('Five_X_final_Late_baseline.mat','x_save')
%%    
%     figure
%     subplot(2,2,1)
%     %plot(1:nGen,n1./(n1+n2+n3+n4+n5))
%     hold on
%     plot(1:nGen,x1)
%     plot(1:nGen,x2)
%     plot(1:nGen,x3)
%     plot(1:nGen,x4)
%     %plot(1:nGen,n2./(n1+n2+n3+n4+n5))
%     %plot(1:nGen,n3./(n1+n2+n3+n4+n5))
%     %plot(1:nGen,n4./(n1+n2+n3+n4+n5))
%     plot(1:nGen,x5)
%     %plot(1:nGen,n5./(n1+n2+n3+n4+n5))
%     title('Late Horizontal')
%     ylim([0 1])
%     legend('1','2','3','4','5')
%     
%     subplot(2,2,2)
%     hold on
%     plot(1:nGen,n1)
%     plot(1:nGen,n2)
%     plot(1:nGen,n3)
%     plot(1:nGen,n4)
%     plot(1:nGen,n5)
%     
%     %title('(2) Early Horizontal')
%     %ylim([0 1])
%     legend('1','2','3','4','5')
%%    
    finalPopSizeH_ph(:) = [n1(1,end); n2(1,end); n3(1,end); n4(1,end); n5(1,end)];
    finalXPropH_ph(:) = [x1(1,end); x2(1,end); x3(1,end); x4(1,end); x5(1,end)];
    finalAgeStructure_ph(:) = [n1(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end)); n2(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end)); n3(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end));n4(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end));n5(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end))];
    
    FREQ_STRAT1(m,:) = finalXPropH_ph;
    
    
    mx1 = mean([x1' x2' x3' x4' x5'],2);
    mn1 = sum([n1' n2' n3' n4' n5'],2);
    %% Early horizontal
    
    %age class frequencies
    n1 = zeros(1,nGen);
    n2 = zeros(1,nGen);
    n3 = zeros(1,nGen);
    n4 = zeros(1,nGen);
    n5 = zeros(1,nGen);
    %informed individuals
    x1 = zeros(1,nGen);
    x2 = zeros(1,nGen);
    x3 = zeros(1,nGen);
    x4 = zeros(1,nGen);
    x5 = zeros(1,nGen);
    propn2 = zeros(1,nGen);
    %first gen
    n1(1) = 0.2*Ni;
    n2(1) = 0.2*Ni;
    n3(1) = 0.2*Ni;
    n4(1) = 0.2*Ni;
    n5(1) = 0.2*Ni;
    x1(1) = 0.001;
    x2(1) = 0.001;
    x3(1) = 0.001;
    x4(1) = 0.001;% trait innovated by reproductive adults
    x5(1) = 0.001;
    
    n1(2) = 0.2*Ni;
    n2(2) = 0.2*Ni;
    n3(2) = 0.2*Ni;
    n4(2) = 0.2*Ni;
    n5(2) = 0.2*Ni;
    x1(2) = 0.001;
    x2(2) = 0.001;
    x3(2) = 0.001;
    x4(2) = 0.001;% trait innovated by reproductive adults
    x5(2) = 0.001;
    
    finalPopSizeH_ph = zeros(5,numel(pv));
    finalXPropH_ph = zeros(5,numel(pv));
    finalAgeStructure_ph = zeros(5,numel(pv));
    
    fiftyPopSizeH_ph = zeros(5,numel(pv));
    fiftyXPropH_ph = zeros(5,numel(pv));
    fiftyAgeStructure_ph = zeros(5,numel(pv));
    
    Nh1000ph = zeros(1,numel(pv));
    
    N = Ni;
    Na=zeros(1,numel(pv));
    %main loop vertical model
    
    %pv = 0.75;%pvi(k);
    for j = 3:nGen
        n1(j) = n4(j-1)*((f4b+w4)*x4(j-1)+f4b*(1-x4(j-1)));
        x1(j) = pv*(n4(j-1)*(f4b+w4)*x4(j-1))/(n4(j-1)*((f4b+w4)*x4(j-1)+f4b*(1-x4(j-1))));%these guys move to age class 2 now
        
        n2(j) = (x1(j-1)*(s1+ws)+(1-x1(j-1))*s1)*n1(j-1);
        x2(j) = x1(j-1)+(1-x1(j-1))*ph*((x1(j-1)*n1(j-1)+x2(j-1)*n2(j-1)+x3(j-1)*n3(j-1)+x4(j-1)*n4(j-1)+x5(j-1)*n5(j-1))/(n1(j-1)+n2(j-1)+n3(j-1)+n4(j-1)+n5(j-1)));
        
        n3(j) = (x2(j-1)*(s2+ws)+(1-x2(j-1))*s2)*n2(j-1);
        x3(j) = x2(j-1)+(1-x2(j-1))*ph*((x2(j-1)*n2(j-1)+x3(j-1)*n3(j-1)+x4(j-1)*n4(j-1)+x5(j-1)*n5(j-1))/(n2(j-1)+n3(j-1)+n4(j-1)+n5(j-1)));
        
        n4(j) = (x3(j-1)*(s3+ws)+(1-x3(j-1))*s3)*n3(j-1);
        x4(j) = x3(j-1)+(1-x3(j-1))*ph*((x3(j-1)*n3(j-1)+x4(j-1)*n4(j-1)+x5(j-1)*n5(j-1))/(n3(j-1)+n4(j-1)+n5(j-1)));
        
        n5(j) = n4(j-1)*x4(j-1)*(s4+ws)+n4(j-1)*(1-x4(j-1))*(s4)+n5(j-1)*x5(j-1)*(s5+ws)+n5(j-1)*(1-x5(j-1))*(s5);
        mixProp = (n4(j-1)*x4(j-1)*(s4+ws)+n4(j-1)*(1-x4(j-1))*(s4))/n5(j);
        x5(j) = x4(j-1)*mixProp+(1-mixProp)*x5(j-1);
        
        N = n1(j)+n2(j)+n3(j)+n4(j)+n5(j);
        
        if N<1
            break
        end
    end
    finalPopSizeH_ph(:) = [n1(1,end); n2(1,end); n3(1,end); n4(1,end); n5(1,end)];
    finalXPropH_ph(:) = [x1(1,end); x2(1,end); x3(1,end); x4(1,end); x5(1,end)];
    finalAgeStructure_ph(:) = [n1(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end)); n2(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end)); n3(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end));n4(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end));n5(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end))];
    
    mx2 = mean([x1' x2' x3' x4' x5'],2);
    mn2 = sum([n1' n2' n3' n4' n5'],2);
    
    FREQ_STRAT2(m,:) = finalXPropH_ph;
    
    %% Hewlett 1 (hg)
    %         %age class frequencies
    %         n1 = zeros(1,nGen);
    %         n2 = zeros(1,nGen);
    %         n3 = zeros(1,nGen);
    %         n4 = zeros(1,nGen);
    %         n5 = zeros(1,nGen);
    %         %informed individuals
    %         x1 = zeros(1,nGen);
    %         x2 = zeros(1,nGen);
    %         x3 = zeros(1,nGen);
    %         x4 = zeros(1,nGen);
    %         x5 = zeros(1,nGen);
    %         %first gen
    %         n1(1) = 0.4*Ni;
    %         n2(1) = 0.2*Ni;
    %         n3(1) = 0.2*Ni;
    %         n4(1) = 0.1*Ni;
    %         n5(1) = 0.1*Ni;
    %         x1(1) = 0.001;
    %         x2(1) = 0.001;
    %         x3(1) = 0.001;
    %         x4(1) = 0.001;% trait innovated by reproductive adults
    %         x5(1) = 0.001;
    %
    %         n1(2) = 0.4*Ni;
    %         n2(2) = 0.2*Ni;
    %         n3(2) = 0.2*Ni;
    %         n4(2) = 0.1*Ni;
    %         n5(2) = 0.1*Ni;
    %         x1(2) = 0.001;
    %         x2(2) = 0.001;
    %         x3(2) = 0.001;
    %         x4(2) = 0.001;% trait innovated by reproductive adults
    %         x5(2) = 0.001;
    %
    %         n1(3) = 0.4*Ni;
    %         n2(3) = 0.2*Ni;
    %         n3(3) = 0.2*Ni;
    %         n4(3) = 0.1*Ni;
    %         n5(3) = 0.1*Ni;
    %         x1(3) = 0.001;
    %         x2(3) = 0.001;
    %         x3(3) = 0.001;
    %         x4(3) = 0.001;% trait innovated by reproductive adults
    %         x5(3) = 0.001;
    %
    %         finalPopSizeH_ph = zeros(5,1);
    %         finalXPropH_ph = zeros(5,1);
    %         finalAgeStructure_ph = zeros(5,1);
    %
    %         Nh1000ph = zeros(1,1);
    %
    %         N = Ni;
    %         Nh=zeros(1,1);
    %         %main loop vertical model
    %
    %
    %         for j = 4:nGen
    %             %n1(j)=n4(j-1)*f4b;
    %             n1(j) = n4(j-1)*((f4b+w4)*x4(j-1)+f4b*(1-x4(j-1)));
    %             x1(j) = pv*n4(j-1)*((f4b+w4)*x4(j-1))/(n4(j-1)*((f4b+w4)*x4(j-1)+f4b*(1-x4(j-1))));%these guys move to age class 2 now
    %
    %             V2 = 0.851;
    %             vert2 = (n4(j-2)*(s4+ws)^2*(f4b+w4)*x4(j-2))/(n4(j-2)*(s4+ws)^2*(f4b+w4)*x4(j-2)+n4(j-2)*f4b*s4^2*(1-x4(j-2)));
    %             n2(j) = (x1(j-1)*(s1+ws)+(1-x1(j-1))*s1)*n1(j-1);
    %             x2(j) = x1(j-1)+(1-x1(j-1))*(V2*vert2*pv+(1-V2)*ph*(((x1(j-1)*n1(j-1)+x2(j-1)*n2(j-1)+x3(j-1)*n3(j-1)+x4(j-1)*n4(j-1)+x5(j-1)*n5(j-1))/(n1(j-1)+n2(j-1)+n3(j-1)+n4(j-1)+n5(j-1)))));
    %
    %             V3 = 0.62;
    %             vert3 = (n4(j-3)*(s4+ws)^3*(f4b+w4)*x4(j-3))/(n4(j-3)*(s4+ws)^3*(f4b+w4)*x4(j-3)+n4(j-3)*f4b*s4^3*(1-x4(j-3)));
    %             n3(j) = (x2(j-1)*(s2+ws)+(1-x2(j-1))*s2)*n2(j-1);%n2(j-1)*s2;
    %             x3(j) = x2(j-1)+(1-x2(j-1))*(V3*vert3*pv+(1-V3)*ph*(((x2(j-1)*n2(j-1)+x3(j-1)*n3(j-1)+x4(j-1)*n4(j-1)+x5(j-1)*n5(j-1))/(n2(j-1)+n3(j-1)+n4(j-1)+n5(j-1)))));
    %
    %             n4(j) = (x3(j-1)*(s3+ws)+(1-x3(j-1))*s3)*n3(j-1);%n3(j-1)*s3;
    %             x4(j) = x3(j-1)+(1-x3(j-1))*ph*((x3(j-1)*n3(j-1)+x4(j-1)*n4(j-1)+x5(j-1)*n5(j-1))/(n3(j-1)+n4(j-1)+n5(j-1)));
    %
    %             n5(j) = n4(j-1)*x4(j-1)*(s4+ws)+n4(j-1)*(1-x4(j-1))*(s4)+n5(j-1)*x5(j-1)*(s5+ws)+n5(j-1)*(1-x5(j-1))*(s5);
    %             mixProp = (n4(j-1)*x4(j-1)*(s4+ws)+n4(j-1)*(1-x4(j-1))*(s4))/n5(j);
    %             x5(j) = x4(j-1)*mixProp+(1-mixProp)*x5(j-1);
    %
    %
    %             N = n1(j)+n2(j)+n3(j)+n4(j)+n5(j);
    %
    %             if n1(1,j)<=1 | n2(1,j)<=1 |n3(1,j)<=1
    %                 break
    %             end
    %         end
    %
    %         finalPopSizeH_ph(:) = [n1(1,end); n2(1,end); n3(1,end); n4(1,end); n5(1,end)];
    %         finalXPropH_ph(:) = [x1(1,end); x2(1,end); x3(1,end); x4(1,end); x5(1,end)];
    %         finalAgeStructure_ph(:) = [n1(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end)); n2(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end)); n3(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end));n4(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end));n5(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end))];
    %
    %         FREQ_STRAT3(m,l,:) = finalXPropH_ph;
    %
    %         mx1 = mean([x1' x2' x3' x4' x5'],2);
    %         mn1 = sum([n1' n2' n3' n4' n5'],2);
    %% Hewlett 2 (ag)
    %
    %         %age class frequencies
    %         n1 = zeros(1,nGen);
    %         n2 = zeros(1,nGen);
    %         n3 = zeros(1,nGen);
    %         n4 = zeros(1,nGen);
    %         n5 = zeros(1,nGen);
    %         %informed individuals
    %         x1 = zeros(1,nGen);
    %         x2 = zeros(1,nGen);
    %         x3 = zeros(1,nGen);
    %         x4 = zeros(1,nGen);
    %         x5 = zeros(1,nGen);
    %         %first gen
    %         n1(1) = 0.4*Ni;
    %         n2(1) = 0.2*Ni;
    %         n3(1) = 0.2*Ni;
    %         n4(1) = 0.1*Ni;
    %         n5(1) = 0.1*Ni;
    %         x1(1) = 0.001;
    %         x2(1) = 0.001;
    %         x3(1) = 0.001;
    %         x4(1) = 0.001;% trait innovated by reproductive adults
    %         x5(1) = 0.001;
    %
    %         n1(2) = 0.4*Ni;
    %         n2(2) = 0.2*Ni;
    %         n3(2) = 0.2*Ni;
    %         n4(2) = 0.1*Ni;
    %         n5(2) = 0.1*Ni;
    %         x1(2) = 0.001;
    %         x2(2) = 0.001;
    %         x3(2) = 0.001;
    %         x4(2) = 0.001;% trait innovated by reproductive adults
    %         x5(2) = 0.001;
    %
    %         n1(3) = 0.4*Ni;
    %         n2(3) = 0.2*Ni;
    %         n3(3) = 0.2*Ni;
    %         n4(3) = 0.1*Ni;
    %         n5(3) = 0.1*Ni;
    %         x1(3) = 0.001;
    %         x2(3) = 0.001;
    %         x3(3) = 0.001;
    %         x4(3) = 0.001;% trait innovated by reproductive adults
    %         x5(3) = 0.001;
    %
    %         finalPopSizeH_ph = zeros(5,1);
    %         finalXPropH_ph = zeros(5,1);
    %         finalAgeStructure_ph = zeros(5,1);
    %
    %         Nh1000ph = zeros(1,1);
    %
    %         N = Ni;
    %         Nh=zeros(1,1);
    %         %main loop vertical model
    %
    %
    %         for j = 4:nGen
    %             %n1(j)=n4(j-1)*f4b;
    %             n1(j) = n4(j-1)*((f4b+w4)*x4(j-1)+f4b*(1-x4(j-1)));
    %             x1(j) = pv*n4(j-1)*((f4b+w4)*x4(j-1))/(n4(j-1)*((f4b+w4)*x4(j-1)+f4b*(1-x4(j-1))));%these guys move to age class 2 now
    %
    %             V2 = 0.48;
    %             vert2 = (n4(j-2)*(s4+ws)^2*(f4b+w4)*x4(j-2))/(n4(j-2)*(s4+ws)^2*(f4b+w4)*x4(j-2)+n4(j-2)*f4b*s4^2*(1-x4(j-2)));
    %             n2(j) = (x1(j-1)*(s1+ws)+(1-x1(j-1))*s1)*n1(j-1);
    %             x2(j) = x1(j-1)+(1-x1(j-1))*(V2*vert2*pv+(1-V2)*ph*(((x1(j-1)*n1(j-1)+x2(j-1)*n2(j-1)+x3(j-1)*n3(j-1)+x4(j-1)*n4(j-1)+x5(j-1)*n5(j-1))/(n1(j-1)+n2(j-2)+n3(j-1)+n4(j-1)+n5(j-1)))));
    %
    %             V3 = 0.083;
    %             vert3 = (n4(j-3)*(s4+ws)^3*(f4b+w4)*x4(j-3))/(n4(j-3)*(s4+ws)^3*(f4b+w4)*x4(j-3)+n4(j-3)*f4b*s4^3*(1-x4(j-3)));
    %             n3(j) = (x2(j-1)*(s2+ws)+(1-x2(j-1))*s2)*n2(j-1);%n2(j-1)*s2;
    %             x3(j) = x2(j-1)+(1-x2(j-1))*(V3*vert3*pv+(1-V3)*ph*(((x2(j-1)*n2(j-1)+x3(j-1)*n3(j-1)+x4(j-1)*n4(j-1)+x5(j-1)*n5(j-1))/(n2(j-1)+n3(j-1)+n4(j-1)+n5(j-1)))));
    %
    %             n4(j) = (x3(j-1)*(s3+ws)+(1-x3(j-1))*s3)*n3(j-1);%n3(j-1)*s3;
    %             x4(j) = x3(j-1)+(1-x3(j-1))*ph*((x3(j-1)*n3(j-1)+x4(j-1)*n4(j-1)+x5(j-1)*n5(j-1))/(n3(j-1)+n4(j-1)+n5(j-1)));
    %
    %             n5(j) = n4(j-1)*x4(j-1)*(s4+ws)+n4(j-1)*(1-x4(j-1))*(s4)+n5(j-1)*x5(j-1)*(s5+ws)+n5(j-1)*(1-x5(j-1))*(s5);
    %             mixProp = (n4(j-1)*x4(j-1)*(s4+ws)+n4(j-1)*(1-x4(j-1))*(s4))/n5(j);
    %             x5(j) = x4(j-1)*mixProp+(1-mixProp)*x5(j-1);
    %
    %
    %             N = n1(j)+n2(j)+n3(j)+n4(j)+n5(j);
    %
    %             if n1(1,j)<=1 | n2(1,j)<=1 |n3(1,j)<=1
    %                 break
    %             end
    %         end
    %
    %         finalPopSizeH_ph(:) = [n1(1,end); n2(1,end); n3(1,end); n4(1,end); n5(1,end)];
    %         finalXPropH_ph(:) = [x1(1,end); x2(1,end); x3(1,end); x4(1,end); x5(1,end)];
    %         finalAgeStructure_ph(:) = [n1(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end)); n2(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end)); n3(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end));n4(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end));n5(end)/(n1(end)+n2(end)+n3(end)+n4(end)+n5(end))];
    %
    %         FREQ_STRAT4(m,l,:) = finalXPropH_ph;
    %
    %         mx1 = mean([x1' x2' x3' x4' x5'],2);
    %         mn1 = sum([n1' n2' n3' n4' n5'],2);
    %%
    
    %x_save = [x1; x2; x3; x4; x5];
    %save('Five_X_final_Early_baseline.mat','x_save')
    
%     subplot(2,2,3)
%     %plot(1:nGen,n1./(n1+n2+n3+n4+n5))
%     hold on
%     plot(1:nGen,x1)
%     plot(1:nGen,x2)
%     plot(1:nGen,x3)
%     plot(1:nGen,x4)
%     %plot(1:nGen,n2./(n1+n2+n3+n4+n5))
%     %plot(1:nGen,n3./(n1+n2+n3+n4+n5))
%     %plot(1:nGen,n4./(n1+n2+n3+n4+n5))
%     plot(1:nGen,x5)
%     %plot(1:nGen,n5./(n1+n2+n3+n4+n5))
%     %title('(2) Early Horizontal')
%     ylim([0 1])
%     legend('1','2','3','4','5')
%     
%     
%     
%     subplot(2,2,4)
%     hold on
%     plot(1:nGen,n1)
%     plot(1:nGen,n2)
%     plot(1:nGen,n3)
%     plot(1:nGen,n4)
%     plot(1:nGen,n5)
%     
%     title('Early Horizontal')%ylim([0 1])
%     legend('1','2','3','4','5')
%     %legend('late horizontal','early horizontal')
%%
end

figure

plot(w4i,mean(FREQ_STRAT1,2)-mean(FREQ_STRAT1(1,:)))
xlabel('wf')
ylabel('freq')
hold on 
plot(w4i,mean(FREQ_STRAT2,2)-mean(FREQ_STRAT2(1,:)))
xlabel('wf')
ylabel('freq')
title('Fertility benefit')

legend('late H','early H')
