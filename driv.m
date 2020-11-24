% basic lorenz coefficients
sigma = 10;
beta = 8/3;
rho = 28;
f = @(t,a) [
    -sigma*a(1) + sigma*a(2);
    rho*a(1) - a(2) - a(1)*a(3);
    -beta*a(3) + a(1)*a(2);
    ];

% variables
IC = [1 1 1];
step = 0.01;
start = 0;
finish = 100;
iterations = (finish - start)/step;
members = 10;
% Change H to observe different relationships
H = [0 0 0; 0 1 0; 0 0 1];
dim = width(H);
observed = height(H);
Gamma = 0.1;
Xi = 0.1;
Eta = 0.1;

% hatsx
vhats = zeros(iterations, members, dim);
mhats = zeros(iterations, dim);
Chats = zeros(iterations, dim, dim);
Sarr = zeros(iterations, observed, observed);
Karr = zeros(iterations, dim, observed);
varr = zeros(iterations, members, dim);
yarr = zeros(iterations, members, observed);

% Raw lorenz63 solutions
[traw, atrueraw] = ode45(f, start:step:finish, IC);
currval = IC;
araw = zeros(iterations, dim);

% RMSE
RMSE = 0;

for val=1:iterations
    currval = onestep(f, step, currval);
    araw(val,:) = currval;
end

% Initialize ensembles
v0 = zeros(members, dim);
for n=1:members
   v0(n, :) = araw(1,:) + mvnrnd(zeros(1,dim),Gamma*ones(1,dim));
end

%plot3(araw(:,1), araw(:,2), araw(:,3));
%plot3(atrueraw(:,1), atrueraw(:,2), atrueraw(:,3))

initcon = v0;
for j=1:iterations
    % Prediction
    for n=1:members
        %[~, tmp] = ode45(f, [j*step j*step+step], initcon(n,:));
        phi = onestep(f, step, initcon(n,:));
        vhats(j, n, :) = phi + mvnrnd(zeros(size(phi)), Xi*ones(1,dim));
    end

    mhats(j, :) = sum(vhats(j,:,:))/members;

    Csum = zeros(dim, dim);
    for n=1:members
        Csum = Csum + (squeeze(vhats(j,n,:)-mhats(j,:)))*(squeeze(vhats(j,n,:)-mhats(j,:)))';
    end
    Chats(j,:,:) = Csum / (members - 1);

    % Analysis
    S(j,:,:) = H * squeeze(Chats(j, :, :)) * H' + Gamma;
    K(j,:,:) = squeeze(Chats(j, :, :)) * H' * pinv(squeeze(S(j,:,:)), 1e-8);
    for n=1:members
        yarr(j,n,:) = (H * araw(j,:)')' + mvnrnd(zeros(1,observed), Eta * ones(1,observed));
        varr(j,n,:) = (eye(dim)-squeeze(K(j,:,:))*H) * squeeze(vhats(j,n,:)) + squeeze(K(j,:,:))*squeeze(yarr(j,n,:));
    end
    initcon = squeeze(varr(j,:,:));
end

%plot3(varr(:,1,1), varr(:,1,2), varr(:,1,3))
plot3(mhats(:,1), mhats(:,2), mhats(:,3))

for n=1:members
    RMSE = RMSE + sqrt(mean((mean(mean(varr)) - mean(mean(vhats))).^2));
end

RMSE = RMSE / members;
