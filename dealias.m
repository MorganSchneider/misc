% Written by Morgan Schneider


function v_new = dealias(v, dat)

v_new = v;

for k = 1:size(v,1)
    
%     fold_inds = [];
%     diffmat = zeros(1,length(v(k,:)));
%     for n = 1:length(v(k,:))-1
%         diff = v(k,n+1) - v(k,n);
%         if abs(diff) >= 24
%             fold_inds = [fold_inds, n, n+1];
%         end
%     end
%     gradv = diffmat;
%     vline = v(k,:);
%     gradv = gradient(vline);
%     fold_inds = find(abs(gradv) >= 10);
%     inds = [1, fold_inds, size(v,2)];
    
    % this needs to be completely different yikes
    
    figure(3)
    plot(1:size(v,2), v(k,:))
    xlabel('n')
    ylabel('v')
    shg
    
    prompt = 'Input number of folds: ';
    n = input(prompt);
    
    if n > 0
        for i = 1:n
            fprintf('Click boundaries of regions to unfold:\n')
            [x, ~] = ginput(2);
            x = round(x);
            if mean(v(k, x(1):x(2))) < 0
                v_new(k, x(1):x(2)) = v(k, x(1):x(2)) + 2 * dat.params.va;
            elseif mean(v(k, x(1):x(2))) > 0
                v_new(k, x(1):x(2)) = v(k, x(1):x(2)) - 2 * dat.params.va;
            end
        end
    end    
end   
    
    
    
    %% ignore
    
    
    
%     if length(fold_inds) == 2 % if there is only one fold
%         if v(k,1) > 0 && v(k,end) > 0 % if both endpoints are red
%             for m = inds(1):inds(2) % from start of line to end of fold
%                 v_new(k,m) = v(k,m) - 2 * dat.params.va; % unfold velocity DOWN
%             end
%         elseif v(k,1) < 0 && v(k,end) < 0 % if both endpoints are green
%             for m = inds(3):inds(end) % from start of fold to end of line
%                 v_new(k,m) = v(k,m) + 2 * dat.params.va; % unfold velocity UP
%             end
%         end
%     else % if there is more than one fold
%         if mod(length(fold_inds),4) ~= 0 % if one endpoint is folded
%             if v(k,1) > 0 % if leftmost point is folded (left is red)
%                 for m = 1:4:length(inds)-1 % first fold starts at index 1 and ends at index 2, 5/6, 9/10, etc
%                     if gradv(inds(m+1)) < 0 % if velocity goes from red to green at end of fold
%                         for n = inds(m):inds(m+1) % from start of fold to end of fold
%                             v_new(k,n) = v(k,n) - 2 * dat.params.va; % unfold velocity DOWN
%                         end
%                     elseif gradv(inds(m+1)) > 0 % if velocity goes from green to red at end of fold
%                         for n = inds(m):inds(m+1) % from start of fold to end of fold
%                             v_new(k,n) = v(k,n) + 2 * dat.params.va; % unfold velocity UP
%                         end
%                     end
%                 end
%             elseif v(k,1) < 0 % if rightmost point is folded (left is green)
%                 for m = 3:4:length(inds)-1 % first fold starts at index 3 and ends at index 4, 7/8, etc until end of row
%                     if gradv(inds(m)) > 0 % if velocity goes from green to red at start of fold
%                         for n = inds(m):inds(m+1) % from start of fold to end of fold
%                             v_new(k,n) = v(k,n) - 2 * dat.params.va; % unfold velocity DOWN
%                         end
%                     elseif gradv(inds(m)) < 0 % if velocity goes from red to green at start of fold
%                         for n = inds(m):inds(m+1) % from start of fold to end of fold
%                             v_new(k,n) = v(k,n) + 2 * dat.params.va; % unfold velocity UP
%                         end
%                     end
%                 end
%             end
%         elseif mod(length(fold_inds),4) == 0 % if neither or both endpoints are folded
%             if v(k,1) > 0 && v(k,end) < 0 % if both endpoints are folded (left is red and right is green)
%                 for m = 1:4:length(inds)-1 %first fold starts at index 1 and ends at index 2, 5/6, 9/10, etc until end of row
%                     if gradv(inds(m+1)) < 0 || gradv(inds(m)) > 0 % if velocity goes from red to green at end of fold, or green to red at start of fold
%                         for n = inds(m):inds(m+1) % from start of fold to end of fold
%                             v_new(k,n) = v(k,n) - 2 * dat.params.va; % unfold velocity DOWN
%                         end
%                     elseif gradv(inds(m+1)) > 0 || gradv(inds(m)) < 0 % if velocity goes from green to red at end of fold, or red to green at start of fold
%                         for n = inds(m):inds(m+1) % from start of fold to end of fold
%                             v_new(k,n) = v(k,n) + 2 * dat.params.va; % unfold velocity UP
%                         end
%                     end
%                 end
%             else % if neither endpoint is folded (left is green and right is red)
%                 for m = 3:4:length(inds)-1 % first fold starts at index 3 and ends at index 4, 7/8, etc
%                     if gradv(inds(m)) > 0 || gradv(inds(m+1)) < 0 % if velocity goes from green to red at start of fold, or red to green at end of fold
%                         for n = inds(m):inds(m+1) % from start of fold to end of fold
%                             v_new(k,n) = v(k,n) - 2 * dat.params.va; % unfold velocity DOWN
%                         end
%                     elseif gradv(inds(m)) < 0 || gradv(inds(m+1)) > 0 % if velocity goes from red to green at start of fold, or green to red at end of fold
%                         for n = inds(m):inds(m+1) % from start of fold to end of fold
%                             v_new(k,n) = v(k,n) + 2 * dat.params.va; % unfold velocity UP
%                         end
%                     end
%                 end
%             end
%         else
%             figure(3)
%             plot(1:size(v,2), v(k,:))
%             xlabel('n')
%             ylabel('v')
%             shg
%             
%             prompt = 'Input number of folds: ';
%             n = input(prompt);
%             
%             if n > 0
%                 for i = 1:n
%                     fprintf('Click boundaries of regions to unfold:\n')
%                     [x, ~] = ginput(2);
%                     x = round(x);
%                     if mean(v(k, x(1):x(2))) < 0
%                         v_new(k, x(1):x(2)) = v(k, x(1):x(2)) + 2 * dat.params.va;
%                     elseif mean(v(k, x(1):x(2))) > 0
%                         v_new(k, x(1):x(2)) = v(k, x(1):x(2)) - 2 * dat.params.va;
%                     end
%                 end
%             end
%         end
%     end
%     
%     
%     % sanity check
% 
%     if ~isempty(find(abs(gradient(v_new(k,:))) >= 10, 1))
%         figure(3)
%         plot(1:size(v_new,2), v_new(k,:))
%         xlabel('n')
%         ylabel('v')
%         shg
%         
%         prompt = 'Input number of folds: ';
%         n = input(prompt);
%         
%         if n > 0
%             for i = 1:n
%                 fprintf('Click boundaries of regions to unfold:\n')
%                 [x, ~] = ginput(2);
%                 x = round(x);
%                 if mean(v(k, x(1):x(2))) < 0
%                     v_new(k, x(1):x(2)) = v_new(k, x(1):x(2)) + 2 * dat.params.va;
%                 elseif mean(v(k, x(1):x(2))) > 0
%                     v_new(k, x(1):x(2)) = v_new(k, x(1):x(2)) - 2 * dat.params.va;
%                 end
%             end
%         end
%     end
% 
% 
%     
%     
%     
%     figure(1)
%     subplot(3,1,1)
%     plot(1:size(v,2), v(k,:))
%     xlabel('index')
%     ylabel('v')
%     title('v old')
%     subplot(3,1,2)
%     plot(1:size(v,2), v_new(k,:))
%     xlabel('index')
%     ylabel('v')
%     title('v new')
%     subplot(3,1,3)
%     plot(1:size(v,2), gradv)
%     xlabel('index')
%     title('gradv')
%     shg
% 
%     
% end

