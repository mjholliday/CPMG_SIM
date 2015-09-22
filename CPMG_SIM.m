# -*-Octave-*-

function [dat,vcpmg,track_t,track] = FGP_cpmg_sim(K,Cyp_Conc=1e-3,FGP_Conc=6e-3,w_bound=[45,90],time_T2=20,ncyc=[2:2:20],SIMPLE_CPMG=0,dead_enzyme=0,Natoms=1e5) 

kab=K(1);
kba=K(2);
kbc=K(3);
kcb=K(4);
kcd=K(5);
kdc=K(6);

%Fid parameters
Nfid=2^10;
tfid=1;
ftime=(0:tfid/(Nfid-1):tfid);
zf=2^14;

N=Natoms; %Number of atoms to simulate
tau=[time_T2./(4.*ncyc),0]; %time btwn cpmg pulses

w=[0,w_bound]./1000; %put w into ms^-1

if(SIMPLE_CPMG)
  kab=kba=kcd=kdc=0;
  dt=1./[1e-8,1e-8;kbc,1e-8;1e-8,kcb].*1000;
  Pop=[0,1/(1+kbc/kcb),1/(1+kcb/kbc)];
  Pop/=sum(Pop);
  EQ_conc=Cyp_Conc.*[0,1/(1+kbc/kcb),1/(1+kcb/kbc),0,0];
else
  %Solve for equilibruim conditions
  ic=[FGP_Conc*0.01,FGP_Conc*.5,FGP_Conc*.5,FGP_Conc*0.01,Cyp_Conc-1e-3];
  [f,EQ_conc,cvg,iter]=leasqr([kab,kba,kbc,kcb,kcd,kdc],[0,0,0,Cyp_Conc,FGP_Conc],ic,"solve_EQ",1e-4,1e6,ones(1,5),0.001*ones(1,5),"partial_EQ");

  kex_un=(kba*kdc*EQ_conc(4)+kcd*kab*EQ_conc(1))/(kdc*EQ_conc(4)+kab*EQ_conc(1));


  if(dead_enzyme) 

    if(dead_enzyme==1)
      per_free=12.7;
    else
      per_free=dead_enzyme;
    end
    
    kad=per_free/(100-per_free)*100;
    kda=1*100;
    
    kbc=kcb=0;
    ic=[FGP_Conc*0.01,FGP_Conc*.5,FGP_Conc*.5,FGP_Conc*0.01,Cyp_Conc-1e-3];
    [f,EQ_conc,cvg,iter]=leasqr([kab,kba,kbc,kcb,kcd,kdc,kad,kda],[0,0,0,Cyp_Conc,FGP_Conc],ic,"solve_EQ_nocat",1e-4,1e6,ones(1,5),0.001*ones(1,5));

  end

   EQ_conc;
  K=[kab*EQ_conc(1),kdc*EQ_conc(4);kbc,kba;kcd,kcb];
  dt=1./K.*1000;
  Pop=[EQ_conc(5),EQ_conc(2),EQ_conc(3)]./(EQ_conc(5)+EQ_conc(2)+EQ_conc(3));
end

  Sout=zeros(size(ncyc,2)+1,1);

%track positon of atom #1
track=[];
track_t=[];

I=zeros(size(ncyc,2)+1,1);

for n=1:size(ncyc,2)+1;
  S=ones(N,1); %Magnitization
  s=ones(N,1); %state (1,2,3)
  t=zeros(N,1); %time
  tc=zeros(N,1)+tau(n); %cpmg time
 
  s(floor(N*Pop(1))+1:N)++;
  s(floor(N*(Pop(1)+Pop(2)))+1:N)++;

   ts = [dt(:,1)(s),dt(:,2)(s)].*rande(N,2);
   sd_both = (ts==repmat(min(ts'),2,1)'); %direction of switch
   sd = sd_both(:,1)*1+sd_both(:,2)*-1; % 1 == forward(ie state 1->2 or 3->1) -1== backwards(ie 3->2 or 1->3)  
   ts = min(ts')'; %switchs in first time point to come up

   if(n==1)
     track=s(1);
     track_t=0;
   end

  if(n<=size(ncyc,2));
    do
      sc=(t+ts<tc); %sc=atoms that will change state before next cpmg pulse
      cp=((!sc+2*(tc>=time_T2))==1); %cp=atoms that will hit cpmg before next state change and are not at end
      
      t(sc)+=ts(sc);
      pro2=((sc+(s==2))==2); %Atoms in state 2 that will switch before next cpmg point
      pro3=((sc+(s==3))==2); %Atoms in state 3 that will switch before next cpmg point
      S(pro2)=S(pro2).*exp(2*pi()*i*ts(pro2)*w(2)); %evolve those in state 2
      S(pro3)=S(pro3).*exp(2*pi()*i*ts(pro3)*w(3)); %evolve those in state 3
     
      s(sc)+=sd(sc); %move to next state
      s(s==0)=3;   %if moved from 1->3
      s(s==4)=1;   %if moved from 3->1
      ts_temp = [dt(:,1)(s(sc)),dt(:,2)(s(sc))].*rande(sum(sc),2); %find new times until switch
      sd_both = (ts_temp==repmat(min(ts_temp'),2,1)'); %direction of switch
      sd(sc) = sd_both(:,1)*1+sd_both(:,2)*-1; % 1 == forward(ie state 1->2 or 3->1) -1== backwards(ie 3->2 or 1->3)  
      ts(sc) = min(ts_temp')'; %switchs in first time point to come up
   
      
      if(n==1)
        if(sc(1)==1)
          if(size(track,1)==0) 
            track=[track,s(1)];
            track_t=[track_t,t(1)];
          else
            track=[track,track(end),s(1)];
            track_t=[track_t,t(1),t(1)];
          end
        end
      end
     
      pro2=((cp+(s==2))==2);  % Atoms in state 2 that reach next cpmg point
      pro3=((cp+(s==3))==2);  % Atoms in state 3 that reach next cpmg point
      S(pro2)=S(pro2).*exp(2*pi()*i*(tc(pro2).-t(pro2)).*w(2)); %evolve state 2
      S(pro3)=S(pro3).*exp(2*pi()*i*(tc(pro3).-t(pro3)).*w(3)); %evolve state 3
      S(cp,:)=conj(S(cp,:));

      ts(cp)-=tc(cp)-t(cp);
      t(cp)=tc(cp);
      tc(cp)+=2*tau(n);
      
      tc(tc>time_T2)=time_T2;
    
    until(!cp && !sc)
    track=[track,track(end)];
    track_t=[track_t,time_T2];
    pro2=(s==2);
    pro3=(s==3);
    S(pro2)=S(pro2).*exp(2*pi()*i*(time_T2-t(pro2))*w(2));
    S(pro3)=S(pro3).*exp(2*pi()*i*(time_T2-t(pro3))*w(3));
    ts-=time_T2-t;
  end
  Sout(n)=sum(S);

%Calculate lineshapes:
%w_fid=[0,w_bound]*10;
%R2free=R2trans=R2cis=40;
%M=[i*2*pi()*w_fid(1)-R2free-kab*EQ_conc(1)-kdc*EQ_conc(4) ,  kba , kcd ;\
%kab*EQ_conc(1) , 2*pi()*i*w_fid(2)-R2trans-kba-kbc , +kcb ;\
%kdc*EQ_conc(4) , kbc , 2*pi()*i*w_fid(3)-R2cis-kcd-kcb];
%
%%Calculate Eigenvalues/vectors
%%rcond(M)
%[V,lam]=eig(M);
%%Set initial conditions/Calculate pre-exponential terms for specific sol'n
%M0=[sum(S(s==1));sum(S(s==2));sum(S(s==3))];
%C=inv(V)*M0;
%fid=sum(V*(exp(diag(lam)*ftime).*repmat(C,1,size(ftime,2))));
%spectrum=real(fft(fid,zf));
%freq=[0:1/t(2)/(size(spectrum,2)):1/t(2)-t(2)/(size(spectrum,2))];

I(n)=abs(Sout(n));%max(spectrum);

end

%I/=max(I);
dat=-(1/(0.001*time_T2))*log(I(1:end-1)/I(end));
vcpmg=1000*(1./(4*tau(1:end-1)));

end
