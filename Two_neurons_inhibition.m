%Constants
Cm=0.01; % Membrane Capacitance
dt=0.04;
t=0:dt:100;
ENa=55.17; % Na reversal potential
EK=-72.14; % K reversal potential
El=-49.42; % Leakage reversal potential
g_Na=1.2; % Na conductance
g_K=0.36; % K conductance
g_l=0.003; % Leakage conductance

A1=.05; %controls the strength of inhibition
A2=.07;

%Initial values neuron 1
V(1)=-60; % Membrane voltage
alpha_m(1)=0.1*(V(1)+35)/(1-exp(-(V(1)+35)/10));
beta_m(1)=4.0*exp(-0.0556*(V(1)+60));
alpha_n(1)=0.01*(V(1)+50)/(1-exp(-(V(1)+50)/10));
beta_n(1)=0.125*exp(-(V(1)+60)/80);
alpha_h(1)=0.07*exp(-0.05*(V(1)+60));
beta_h(1)=1/(1+exp(-(0.1)*(V(1)+30)));

stimulus1(1)=.7; %External current
I1(1)=0; %inhibition
time1=0; %time since neuron 1 has fired

m(1)=alpha_m(1)/(alpha_m(1)+beta_m(1));
n(1)=alpha_n(1)/(alpha_n(1)+beta_n(1));
h(1)=alpha_h(1)/(alpha_h(1)+beta_h(1));


%Initial values neuron 2
V2(1)=-60; % Membrane voltage
alpha2_m(1)=0.1*(V2(1)+35)/(1-exp(-(V2(1)+35)/10));
beta2_m(1)=4.0*exp(-0.0556*(V2(1)+60));
alpha2_n(1)=0.01*(V2(1)+50)/(1-exp(-(V2(1)+50)/10));
beta2_n(1)=0.125*exp(-(V2(1)+60)/80);
alpha2_h(1)=0.07*exp(-0.05*(V2(1)+60));
beta2_h(1)=1/(1+exp(-(0.1)*(V2(1)+30)));

stimulus2(1)=.5;
I2(1)=0;
time2=0; %time since neuron 2 has fired

m2(1)=alpha2_m(1)/(alpha2_m(1)+beta2_m(1));
n2(1)=alpha2_n(1)/(alpha2_n(1)+beta2_n(1));
h2(1)=alpha2_h(1)/(alpha2_h(1)+beta2_h(1));

%begin loop
for i=1:length(t)-1

 %Euler method for m,n,h
 m(i+1)=m(i)+dt*((alpha_m(i)*(1-m(i)))-(beta_m(i)*m(i)));
 n(i+1)=n(i)+dt*((alpha_n(i)*(1-n(i)))-(beta_n(i)*n(i)));
 h(i+1)=h(i)+dt*((alpha_h(i)*(1-h(i)))-(beta_h(i)*h(i)));
 
 m2(i+1)=m2(i)+dt*((alpha2_m(i)*(1-m2(i)))-(beta2_m(i)*m2(i)));
 n2(i+1)=n2(i)+dt*((alpha2_n(i)*(1-n2(i)))-(beta2_n(i)*n2(i)));
 h2(i+1)=h2(i)+dt*((alpha2_h(i)*(1-h2(i)))-(beta2_h(i)*h2(i)));

 %Euler method for V
 INa=g_Na*m(i)^3*h(i)*(V(i)-ENa);
 IK=g_K*n(i)^4*(V(i)-EK);
 Il=g_l*(V(i)-El);
 
 INa2=g_Na*m2(i)^3*h2(i)*(V2(i)-ENa);
 IK2=g_K*n2(i)^4*(V2(i)-EK);
 Il2=g_l*(V2(i)-El);
 
 V(i+1)=V(i)+(dt)*((1/Cm)*(I1(1)+stimulus1(i)-(INa+IK+Il)));
 
 V2(i+1)=V2(i)+(dt)*((1/Cm)*(I2(1)+stimulus2(i)-(INa2+IK2+Il2)));
 
 %check if either neuron fires
 if V(i+1)>0 && V(i)<=0
    I2(i+1)=I2(i)-A2; %if neuron 1 has fired, inhibit neuron 2
    time1=0; 
 else
    time1=time1+dt;
    I2(i+1)=I2(i)*exp(-time1/5); %otherwise, exponential decay of the inhibition term
 end
  
 if V2(i+1)>0 && V2(i)<=0
    I1(i+1)=I1(i)-A1;
    time2=0;
 else
    time2=time2+dt;
    I1(i+1)=I1(i)*exp(-time2/5);
 end
 
 %recalculating alpha and beta values
 alpha_m(i+1)=0.1*(V(i+1)+35)/(1-exp(-(V(i+1)+35)/10));
 beta_m(i+1)=4.0*exp(-0.0556*(V(i+1)+60));
 alpha_n(i+1)=0.01*(V(i+1)+50)/(1-exp(-(V(i+1)+50)/10));
 beta_n(i+1)=0.125*exp(-(V(i+1)+60)/80);
 alpha_h(i+1)=0.07*exp(-0.05*(V(i+1)+60));
 beta_h(i+1)=1/(1+exp(-(0.1)*(V(i+1)+30)));
 
 alpha2_m(i+1)=0.1*(V2(i+1)+35)/(1-exp(-(V2(i+1)+35)/10));
 beta2_m(i+1)=4.0*exp(-0.0556*(V2(i+1)+60));
 alpha2_n(i+1)=0.01*(V2(i+1)+50)/(1-exp(-(V2(i+1)+50)/10));
 beta2_n(i+1)=0.125*exp(-(V2(i+1)+60)/80);
 alpha2_h(i+1)=0.07*exp(-0.05*(V2(i+1)+60));
 beta2_h(i+1)=1/(1+exp(-(0.1)*(V2(i+1)+30)));
 
 %update stimulus (white noise) values
 stimulus1(i+1)=normrnd(.5,.6);
 stimulus2(i+1)=normrnd(.5,.6);
 
end

figure(1)
plot(t,V,t,V2)
xlabel('Time (ms)')
ylabel('Voltage (mv)')
title('Voltage over time')
legend('Neuron 1','Neuron 2')

figure(2)
plot(t,I1,t,I2)
title('Inhibition over time')

figure(3)
plot(t,m,t,n,t,h)
title('Gating variables over time')
legend('m','n','h')