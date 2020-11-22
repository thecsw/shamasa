% basic lorenz coefficients
sigma = 10;
beta = 8/3;
rho = 28;
f = @(t,a) [
    -sigma*a(1) + sigma*a(2);
    rho*a(1) - a(2) - a(1)*a(3);
    -beta*a(3) + a(1)*a(2)];

% variables
IC = [1 1 1];
start = 0;
finish = 100;
iterations = 42;
members = 10;
step = 0.01;
% Change H to observe different relationships
H = [1 0 0; 0 1 0; 0 0 1];
dim = height(H);
observed = width(H);
Gamma = 0.000001;
Xi = 0.0001;
Eta = 0.0001;

% hats
vhats = zeros(iterations, members, dim);
mhats = zeros(iterations, dim);
Chats = zeros(iterations, dim, dim);
Sarr = zeros(iterations, observed, observed);
Karr = zeros(iterations, dim, observed);
varr = zeros(iterations, members, dim);
yarr = zeros(iterations, members, observed);

% Raw lorenz63 solutions
[traw, araw] = ode45(f, start:step:finish, IC);

% Initialize ensembles
v0 = zeros(members, dim);
for n=1:members
    v0(n, :) = araw(1,:) + mvnrnd(zeros(1,3),Gamma*ones(1,3));
end

initcon = v0;
for j=1:iterations
    % Prediction
    for n=1:members
        [~, tmp] = ode45(f, [0 step/2 step], initcon(n,:));

        vhats(j, n, :) = tmp(end,:);
        %+ mvnrnd(zeros(size(tmp(end,:))), Xi*ones(1,dim));
    end
    
    mhats(j, :) = sum(vhats(j,:,:))/members;
    
    Csum = zeros(dim, dim);
    for n=1:members
        Csum = Csum + (squeeze(vhats(j,n,:)-mhats(j,:)))*(squeeze(vhats(j,n,:)-mhats(j,:)))';
    end
    Chats(j,:,:) = Csum / (members - 1);
    
    % Analysis
    S(j,:,:) = H * squeeze(Chats(j, :, :)) * H' + Gamma;
    K(j,:,:) = squeeze(Chats(j, :, :)) * H' * inv(squeeze(S(j,:,:)));
    for n=1:members
        yarr(j,n,:) = (H * araw(j,:)')' + mvnrnd(zeros(1,observed), Eta * ones(1,observed));
        varr(j,n,:) = (eye(dim)-squeeze(K(j,:,:))*H) * squeeze(vhats(j,n,:)) + squeeze(K(j,:,:))*squeeze(yarr(j,n,:));
    end

    initcon = squeeze(varr(j,:,:));
end

%plot3(araw(:,1), araw(:,2), araw(:,3))
%plot3(varr(:,1,1), varr(:,1,2), varr(:,1,3))
plot3(mhats(:,1), mhats(:,2), mhats(:,3))


