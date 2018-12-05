% 64BIT 4个1的熵与期望计算

a = [560,25088,21952,301056,286720];
b = [7,0,0,0,1;6,1,0,1,0;6,0,2,0,0;5,2,1,0,0;4,4,0,0,0];

cnt(1,1) = sum(a.*b(:,1)');
cnt(1,2) = sum(a.*b(:,2)');
cnt(1,3) = sum(a.*b(:,3)');
cnt(1,4) = sum(a.*b(:,4)');
cnt(1,5) = sum(a.*b(:,5)');
cnt(2,:) = cnt(1,:)/sum(cnt(1,:));
HX = 0;
totalcnt = sum(cnt(1,:));
for i = 1:length(cnt)
    HX = HX-cnt(2,i)*cnt(1,i)/(totalcnt*cnt(2,i))*log2(cnt(1,i)/(totalcnt*cnt(2,i)));
end
cnt(3,:) = cnt(1,:)/totalcnt./cnt(2,:);
alphabet = [];
for i = 1:sum(cnt(2,:))
    str = ['a' num2str(i)];
    alphabet = [alphabet,{str}];
end
prob = [];
for i = 1:length(cnt)
    for j = 1:cnt(2,i)
        prob = [prob cnt(3,i)];
    end
end
   
dict     = huffmandict_(alphabet,prob, 0); % Set 0->1, creates log file.
for i = 1:length(dict.code)
    if str2num(cell2mat(dict.code(1,i))) == 0
        dict.codelength(1,i) = 1;
    else
        dict.codelength(1,i) = length(cell2mat(dict.code(1,i)));
    end
end
EX = sum(dict.codelength.*prob);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%1023bit 9个1的H和E计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入方案设计
Type = [0,0,0,0,0,0,0,0,1,32;1,0,0,0,0,0,0,1,0,0;0,1,0,0,0,0,1,0,0,0;2,0,0,0,0,0,1,0,0,0;1,0,1,0,0,1,0,0,0,0;1,1,0,0,0,1,0,0,0,0;3,0,0,0,0,1,0,0,0,0;0,0,0,1,1,0,0,0,0,0;1,0,1,0,1,0,0,0,0,0;0,2,0,0,1,0,0,0,0,0;2,1,0,0,1,0,0,0,0,0;4,0,0,0,1,0,0,0,0,0;1,0,0,2,0,0,0,0,0,0;0,1,1,1,0,0,0,0,0,0;2,0,1,1,0,0,0,0,0,0;1,2,0,1,0,0,0,0,0,0;3,1,0,1,0,0,0,0,0,0;5,0,0,1,0,0,0,0,0,0;0,0,3,0,0,0,0,0,0,0;1,1,2,0,0,0,0,0,0,0;3,0,2,0,0,0,0,0,0,0;0,3,1,0,0,0,0,0,0,0;2,2,1,0,0,0,0,0,0,0;4,1,1,0,0,0,0,0,0,0;6,0,1,0,0,0,0,0,0,0;1,4,0,0,0,0,0,0,0,0;3,3,0,0,0,0,0,0,0,0;5,2,0,0,0,0,0,0,0,0;7,1,0,0,0,0,0,0,0,0;9,0,0,0,0,0,0,0,0,0;1,2,3,4,5,6,7,8,9,0];

% 各个组对应情况数
% 比如31bit 有2个1，那就有combntns(31,2)个情况
for i = 1:9
    Type(32,i) = combntns(31,i);
end
Type(32,10) = 1;

% 把31bit全部是0的数量算出来
Type(1:30,10) = ones(30,1)*33-sum(Type(1:30,1:9)')';
Type(1:30,11) = 1;

% 对应不同的方案（不同种类组最后组成1023bit 有9个1的信息）各种组又有多少种情况
for i = 1:10
    Type(1:30,11) = Type(1:30,11) .* (Type(32,i).^Type(1:30,i));
end

% 对每个方案组类型确定，各类型数量确定的情况下，计算有多少种排列组合的可能
for i = 1:30
    total = 33-Type(i,10);
    tmp = combntns(33,total);
    for j = 1:9
        tmp = tmp * combntns(total,Type(i,j));
        total = total - Type(i,j);
    end
    Type(i,12) = tmp;
end

% 计算这种方案总数
Type(:,13) = Type(:,12).*Type(:,11);

% 这个是指总共有多少种1023bit的信息
total_type = sum(Type(:,13));
% 基于所有情况，各个类型的组的数量
for i = 1:10
    Type(33,i) = sum(Type(1:30,i).*Type(1:30,13));
end
% 由于比如31bit 1个1有31中情况，Type(33,i)算的是总数，要除以类型数
Type(34,1:10) = Type(33,1:10)./Type(32,1:10);
% 熵计算
HX = 0;
% 这个是指所有情况下，各个31bit分段出现的总数，应为total_type的33倍(一段有33组)
total_typenum = sum(Type(33,1:10));
total_typenum/total_type %验证
for i = 1:10
    HX = HX-Type(32,i)*Type(34,i)/(total_typenum)*log2(Type(34,i)/(total_typenum));
end
% 算得HX = 2.2538
% 中断在这里了，总情况数太多，电脑跑不动
alphabet = [];
for i = 1:sum(Type(32,1:10))
    str = ['a' num2str(i)];
    alphabet = [alphabet,{str}];
end
prob = [];
for i = 1:length(cnt)
    for j = 1:cnt(2,i)
        prob = [prob cnt(3,i)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%511bit 4个1的H和E计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入方案设计
numall = 512;
bitnum = 64;
groupnum = 8;
num1 = 4;
Type = zeros(5,5);
Type = [0,0,0,1,0;1,0,1,0,0;0,2,0,0,0;2,1,0,0,0;4,0,0,0,0];
[row_Type,col_Type] = size(Type);
% 各个组对应情况数
% 比如31bit 有2个1，那就有combntns(31,2)个情况
for i = 1:col_Type-1
    Type(row_Type+1,i) = combntns(bitnum,i);
end
Type(row_Type+1,col_Type) = 1;

% 把组中全部是0的数量算出来
Type(1:row_Type,col_Type) = ones(row_Type,1)*num1-sum(Type(1:row_Type,1:(col_Type-1))')';


% 对应不同的方案（不同种类组最后组成1023bit 有9个1的信息）各种组又有多少种情况
Type(1:row_Type,col_Type+1) = 1;
for i = 1:row_Type
    Type(1:row_Type,col_Type+1) = Type(1:row_Type,col_Type+1) .* (Type(row_Type+1,i).^Type(1:row_Type,i));
end

% 对每个方案组类型确定，各类型数量确定的情况下，计算有多少种排列组合的可能
for i = 1:row_Type
    total = groupnum-Type(i,col_Type);
    tmp = combntns(groupnum,total);
    for j = 1:col_Type-1
        tmp = tmp * combntns(total,Type(i,j));
        total = total - Type(i,j);
    end
    Type(i,col_Type+2) = tmp;
end

% 计算这种方案总数
Type(:,col_Type+3) = Type(:,col_Type+2).*Type(:,col_Type+1);

% 这个是指总共有多少种1023bit的信息
total_type = sum(Type(:,col_Type+3));
% 基于所有情况，各个类型的组的数量
for i = 1:col_Type
    Type(row_Type+2,i) = sum(Type(1:row_Type,i).*Type(1:row_Type,col_Type+3));
end
% 由于比如31bit 1个1有31中情况，Type(33,i)算的是总数，要除以类型数
Type(row_Type+3,1:col_Type) = Type(row_Type+2,1:col_Type)./Type(row_Type+1,1:col_Type);
% 熵计算
HX = 0;
% 这个是指所有情况下，各个31bit分段出现的总数，应为total_type的33倍(一段有33组)
total_typenum = sum(Type(row_Type+2,1:col_Type));
total_typenum/total_type %验证
for i = 1:col_Type
    HX = HX-Type(row_Type+1,i)*Type(row_Type+3,i)/(total_typenum)*log2(Type(row_Type+3,i)/(total_typenum));
end
% 算得HX = 3.8754
alphabet = [];
for i = 1:sum(Type(row_Type+1,1:col_Type))
    str = ['a' num2str(i)];
    alphabet = [alphabet,{str}];
end
prob = [];
for i = 1:col_Type
    for j = 1:Type(row_Type+1,i)
        prob = [prob Type(row_Type+3,i)/total_typenum];
    end
end

dict     = huffmandict_(alphabet,prob, 0); % Set 0->1, creates log file.
for i = 1:length(dict.code)
%     if str2num(cell2mat(dict.code(1,i))) == 0
%         dict.codelength(1,i) = 1;
%     else
        dict.codelength(1,i) = length(cell2mat(dict.code(1,i)));
%     end
end
EX = sum(dict.codelength.*prob);
% 算得EX = 3.9273
