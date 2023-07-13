close all
clear all 
        
nPatches = [20 30 40 50 60 70];% 30 40 50 60 70];%40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200];% 50 60 70 80 90 100];% 110 120 130 140 150 160 170 180 190 200];% 110 120 130 140];
nLoc = [56];% 286];% 364];% 455];
overlap = 0.05;
oversamp = 3;
ep = 0.1;
nu = 0.3;
func = 'gauss';

geom = 'cuboid';

switch geom
    case 'cuboid'
        fldr = '../src/experiments/cuboid/';
    case 'diaphragm'
        fldr = '../src/experiments/diaphragm/';
end    
    

for i = 1:length(nPatches)
    for j = 1:length(nLoc)
        exper_fldr = [fldr, '3D_cuboid_nLoc',int2str(nLoc(j)),'_np0', ...
                    num2str(nPatches(1)),'_Patches', ...
                    int2str(nPatches(i)),'_delta', num2str(overlap),'_q', ...
                    num2str(oversamp),'_ep', num2str(ep),'_nu',int2str(nu*100), ...
                    '_manufactured_',func];

        % pos = find(exper_fldr~='.');
        % exper_fldr = exper_fldr(pos);
        load([exper_fldr,'.mat']);
        L2_RBF(i,j) = Result.L2;
        Linf_RBF(i,j) = Result.Linf;
        db(i,j) = Result.differror.b;
        dxx(i,j) = Result.differror.dxx;
        dxy(i,j) = Result.differror.dxy;
        dxb(i,j) = Result.differror.dxb;
        dzb(i,j) = Result.differror.dzb;        
        [~,dist] = knnsearch(Result.xc,Result.xc,'k',2);
        h(i,j) = max(dist(:,2));
    end
end

t = datetime;
t.Format = 'yyyyMMddHHmmss';

figure()
if length(nPatches) == 1
    [hTemp, id] = sort(h, 'descend');
    loglog(hTemp,dxx(id),'o-');
    hold on
else
    for i = 1:length(nLoc)
        [hTemp, id] = sort(h(:,i), 'descend');
        loglog(hTemp,dxx(id,i),'o-');
        hold on
    end
end
title('dxx')
xlabel('h')
ylabel('error')
savefig(strcat(fldr, string(t), '_differror_dxx'))

figure()
if length(nPatches) == 1
    [hTemp, id] = sort(h, 'descend');
    loglog(hTemp,dxy(id),'o-');
    hold on
else
    for i = 1:length(nLoc)
        [hTemp, id] = sort(h(:,i), 'descend');
        loglog(hTemp,dxy(id,i),'o-');
        hold on
    end
end
title('dxy')
xlabel('h')
ylabel('error')
savefig(strcat(fldr, string(t), '_differror_dxy'))

figure()
if length(nPatches) == 1
    [hTemp, id] = sort(h, 'descend');
    loglog(hTemp,dxb(id),'o-');
    hold on
else
    for i = 1:length(nLoc)
        [hTemp, id] = sort(h(:,i), 'descend');
        loglog(hTemp,dxb(id,i),'o-');
        hold on
    end
end
title('dxb')
xlabel('h')
ylabel('error')
savefig(strcat(fldr, string(t), '_differror_dxb'))

figure()
if length(nPatches) == 1
    [hTemp, id] = sort(h, 'descend');
    loglog(hTemp,db(id),'o-');
    hold on
else
    for i = 1:length(nLoc)
        [hTemp, id] = sort(h(:,i), 'descend');
        loglog(hTemp,db(id,i),'o-');
        hold on
    end
end
title('db')
xlabel('h')
ylabel('error')
savefig(strcat(fldr, string(t), '_differror_db'))

figure()
if length(nPatches) > 1
    for i = 1:length(nLoc)
        [hTemp, id] = sort(h(:,i), 'descend');
        loglog(hTemp,L2_RBF(id,i),'-o','LineWidth',1.5); 
        lgnd{i} = int2str(i+4);
        hold on
    end
else
    [hTemp, id] = sort(h(1,:), 'descend');
    loglog(hTemp,L2_RBF(1,id),'-o','LineWidth',1.5);
    hold on
    lgnd = {'RBF-PUM'};
end

grid on
set(gca,'FontSize',28)
xlabel('h','FontSize',35);
ylabel('||e||_{l_2}','FontSize',35);
hold off
legend(lgnd,'Location','northwest','Orientation','horizontal','FontSize',28)
savefig(strcat(fldr, string(t), '_L2'))


figure()
if length(nPatches) > 1
    for i = 1:length(nLoc)
        [hTemp, id] = sort(h(:,i), 'descend');
        loglog(hTemp,Linf_RBF(id,i),'-o','LineWidth',1.5); 
        lgnd{i} = int2str(i+4);
        hold on
    end
else
    [hTemp, id] = sort(h(1,:), 'descend');
    loglog(hTemp,Linf_RBF(1,id),'-o','LineWidth',1.5);
    hold on
    lgnd = {'RBF-PUM'};
end

grid on
set(gca,'FontSize',28)
xlabel('h','FontSize',35);
ylabel('||e||_{l_inf}','FontSize',35);
hold off
legend(lgnd,'Location','northwest','Orientation','horizontal','FontSize',35)
savefig(strcat(fldr, string(t), '_Linf'))
