function outputs = screation(a)

listofchars = ['.', 'o', 'x', '+', '*', 's', 's', 'd', 'v', '^', '<', '>', 'p', 'h'];
listofcolors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'];
randvals1 = rand(1,length(listofchars));
outputs1 = '.';
% for i = 1:length(listofchars)
%     if randvals1(i)>0.75
%         outputs1 = listofchars(i);
%         break;
%     end
% end


randvals2 = rand(1,length(listofcolors));
outputs2 = 'r';
% for i = 1:length(listofcolors)
%     if randvals2(i)>0.75
%         outputs2 = listofcolors(i);
%         break;
%     end
% end

outputs = strcat(outputs2, outputs1);
