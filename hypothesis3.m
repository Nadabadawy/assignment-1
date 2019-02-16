clear 
ds = datastore('house_data_complete.csv','TreatAsMissing','NA',.....
     'MissingValue',0,'ReadSize',21613);
T = read(ds);
size(T);
Alpha=0.01; % the perfect alpha
m=length(T{:,1});

% taking different features with different polynimials.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross validation 60%
m_60 = ceil(length(T{:,1})*0.6);
m_20 = ceil(length(T{:,1})*0.2);
U0=T{1:m_60,2};
U=T{1:m_60,4:10};
Ux = T{1:m_60,11:19};
U1=T{1:m_60,20:21};
U2=U.^2;
X=[ones(m_60,1) U Ux.^2 U1.^2 Ux.^3];

n=length(X(1,:));
for w=2:n
    if max(abs(X(:,w)))~=0
    X(:,w)=(X(:,w)-mean((X(:,w))))./std(X(:,w));
    end
end

Y=T{1:m_60,3}/mean(T{1:m_60,3});
 Theta=zeros(n,1);
k=1;
lamda=0;

E(k)=(1/(2*m_60))*sum((X*Theta-Y).^2)+(lamda/(2*m_60))*sum((Theta).^2);

R=1;
while R==1
Alpha=Alpha*1;
Theta=Theta-(Alpha/m_60)*X'*(X*Theta-Y);
k=k+1;

E(k)=(1/(2*m_60))*sum((X*Theta-Y).^2);
if E(k-1)-E(k)<0
    break
end 
q=(E(k-1)-E(k))./E(k-1);
if q <.0001;
    R=0;
end
end
figure(1)
plot(E)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross validation 20%
UU0=T{m_60+1:m_60+m_20,2};
UU=T{m_60+1:m_60+m_20,4:10};
UU1=T{m_60+1:m_60+m_20,20:21};
UU2=UU.^2;
UUx = T{m_60+1:m_60+m_20,11:19};
XX=[ones(m_20,1) UU UUx.^2 UU1.^2 UUx.^3];
Theta1=Theta;
YY=T{m_60+1:m_60+m_20,3}/mean(T{m_60+1:m_60+m_20,3});
nn=length(XX(1,:));
kk=1;


for w=2:nn
    if max(abs(XX(:,w)))~=0
    XX(:,w)=(XX(:,w)-mean((XX(:,w))))./std(XX(:,w));
    end
end
EE=(1/(2*m_20))*sum((XX*Theta1-YY).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% last cross validation 20%
p=m-(m_60+m_20);
UUU0=T{m_60+m_20+1:end,2};
UUU=T{m_60+m_20+1:end,4:10};
UUUx=T{m_60+m_20+1:end,11:19};
UUU1=T{m_60+m_20+1:end,20:21};
UUU2=UUU.^2;
XXX=[ones(p,1) UUU UUUx.^2 UUU1.^2 UUUx.^3];
Theta2=Theta1;
YYY=T{m_60+1+m_20:end,3}/mean(T{m_60+1+m_20:end,3});
nnn=length(XXX(1,:));
kkk=1;


for w=2:nnn
    if max(abs(XXX(:,w)))~=0
    XXX(:,w)=(XXX(:,w)-mean((XXX(:,w))))./std(XXX(:,w));
    end
end
EEE(kkk)=(1/(2*p))*sum((XXX*Theta2-YYY).^2);