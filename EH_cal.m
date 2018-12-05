% 64BIT 4��1��������������

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
%%%%%%%%%%%%%%%%%%%%%%%%%1023bit 9��1��H��E����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���뷽�����
Type = [0,0,0,0,0,0,0,0,1,32;1,0,0,0,0,0,0,1,0,0;0,1,0,0,0,0,1,0,0,0;2,0,0,0,0,0,1,0,0,0;1,0,1,0,0,1,0,0,0,0;1,1,0,0,0,1,0,0,0,0;3,0,0,0,0,1,0,0,0,0;0,0,0,1,1,0,0,0,0,0;1,0,1,0,1,0,0,0,0,0;0,2,0,0,1,0,0,0,0,0;2,1,0,0,1,0,0,0,0,0;4,0,0,0,1,0,0,0,0,0;1,0,0,2,0,0,0,0,0,0;0,1,1,1,0,0,0,0,0,0;2,0,1,1,0,0,0,0,0,0;1,2,0,1,0,0,0,0,0,0;3,1,0,1,0,0,0,0,0,0;5,0,0,1,0,0,0,0,0,0;0,0,3,0,0,0,0,0,0,0;1,1,2,0,0,0,0,0,0,0;3,0,2,0,0,0,0,0,0,0;0,3,1,0,0,0,0,0,0,0;2,2,1,0,0,0,0,0,0,0;4,1,1,0,0,0,0,0,0,0;6,0,1,0,0,0,0,0,0,0;1,4,0,0,0,0,0,0,0,0;3,3,0,0,0,0,0,0,0,0;5,2,0,0,0,0,0,0,0,0;7,1,0,0,0,0,0,0,0,0;9,0,0,0,0,0,0,0,0,0;1,2,3,4,5,6,7,8,9,0];

% �������Ӧ�����
% ����31bit ��2��1���Ǿ���combntns(31,2)�����
for i = 1:9
    Type(32,i) = combntns(31,i);
end
Type(32,10) = 1;

% ��31bitȫ����0�����������
Type(1:30,10) = ones(30,1)*33-sum(Type(1:30,1:9)')';
Type(1:30,11) = 1;

% ��Ӧ��ͬ�ķ�������ͬ������������1023bit ��9��1����Ϣ�����������ж��������
for i = 1:10
    Type(1:30,11) = Type(1:30,11) .* (Type(32,i).^Type(1:30,i));
end

% ��ÿ������������ȷ��������������ȷ��������£������ж�����������ϵĿ���
for i = 1:30
    total = 33-Type(i,10);
    tmp = combntns(33,total);
    for j = 1:9
        tmp = tmp * combntns(total,Type(i,j));
        total = total - Type(i,j);
    end
    Type(i,12) = tmp;
end

% �������ַ�������
Type(:,13) = Type(:,12).*Type(:,11);

% �����ָ�ܹ��ж�����1023bit����Ϣ
total_type = sum(Type(:,13));
% ��������������������͵��������
for i = 1:10
    Type(33,i) = sum(Type(1:30,i).*Type(1:30,13));
end
% ���ڱ���31bit 1��1��31�������Type(33,i)�����������Ҫ����������
Type(34,1:10) = Type(33,1:10)./Type(32,1:10);
% �ؼ���
HX = 0;
% �����ָ��������£�����31bit�ֶγ��ֵ�������ӦΪtotal_type��33��(һ����33��)
total_typenum = sum(Type(33,1:10));
total_typenum/total_type %��֤
for i = 1:10
    HX = HX-Type(32,i)*Type(34,i)/(total_typenum)*log2(Type(34,i)/(total_typenum));
end
% ���HX = 2.2538
% �ж��������ˣ��������̫�࣬�����ܲ���
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
%%%%%%%%%%%%%%%%%%%%%%%%%511bit 4��1��H��E����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���뷽�����
numall = 512;
bitnum = 64;
groupnum = 8;
num1 = 4;
Type = zeros(5,5);
Type = [0,0,0,1,0;1,0,1,0,0;0,2,0,0,0;2,1,0,0,0;4,0,0,0,0];
[row_Type,col_Type] = size(Type);
% �������Ӧ�����
% ����31bit ��2��1���Ǿ���combntns(31,2)�����
for i = 1:col_Type-1
    Type(row_Type+1,i) = combntns(bitnum,i);
end
Type(row_Type+1,col_Type) = 1;

% ������ȫ����0�����������
Type(1:row_Type,col_Type) = ones(row_Type,1)*num1-sum(Type(1:row_Type,1:(col_Type-1))')';


% ��Ӧ��ͬ�ķ�������ͬ������������1023bit ��9��1����Ϣ�����������ж��������
Type(1:row_Type,col_Type+1) = 1;
for i = 1:row_Type
    Type(1:row_Type,col_Type+1) = Type(1:row_Type,col_Type+1) .* (Type(row_Type+1,i).^Type(1:row_Type,i));
end

% ��ÿ������������ȷ��������������ȷ��������£������ж�����������ϵĿ���
for i = 1:row_Type
    total = groupnum-Type(i,col_Type);
    tmp = combntns(groupnum,total);
    for j = 1:col_Type-1
        tmp = tmp * combntns(total,Type(i,j));
        total = total - Type(i,j);
    end
    Type(i,col_Type+2) = tmp;
end

% �������ַ�������
Type(:,col_Type+3) = Type(:,col_Type+2).*Type(:,col_Type+1);

% �����ָ�ܹ��ж�����1023bit����Ϣ
total_type = sum(Type(:,col_Type+3));
% ��������������������͵��������
for i = 1:col_Type
    Type(row_Type+2,i) = sum(Type(1:row_Type,i).*Type(1:row_Type,col_Type+3));
end
% ���ڱ���31bit 1��1��31�������Type(33,i)�����������Ҫ����������
Type(row_Type+3,1:col_Type) = Type(row_Type+2,1:col_Type)./Type(row_Type+1,1:col_Type);
% �ؼ���
HX = 0;
% �����ָ��������£�����31bit�ֶγ��ֵ�������ӦΪtotal_type��33��(һ����33��)
total_typenum = sum(Type(row_Type+2,1:col_Type));
total_typenum/total_type %��֤
for i = 1:col_Type
    HX = HX-Type(row_Type+1,i)*Type(row_Type+3,i)/(total_typenum)*log2(Type(row_Type+3,i)/(total_typenum));
end
% ���HX = 3.8754
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
% ���EX = 3.9273
