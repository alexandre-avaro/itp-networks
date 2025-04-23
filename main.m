clear all; close all; clc;

Network_V = 12; % for fixed voltage - not used (set mode in Solve_Network)
Network_I = 0.5e-3;

% Choose Network Properties
Split_Style = 2; % split 1 to Split_Style at each bifurcation
Generation_Count = 5; % split Generation_Count times/generations

generation_lengths = 1e-3 *        [10 7.531 5.302 5.081 5.02 11.505 11.505 5.02 5.081 5.302 7.531 10];
generation_areas = 200e-6 * 1e-6 * [4000 2000 1000 500 250 125 125 250 500 1000 2000 4000]; 
%                         for gens of 1   2    4    8  16  32  32  16   8    4    2    1

% Calculating node & edge count
% add an intermediate set of m^n (i.e. nodes in the middle branches to help with visualization
Num_Edges = 2*(Split_Style^(Generation_Count+1) - 1)/(Split_Style-1);
Num_Nodes = 2*(1 + (Split_Style^Generation_Count - 1)/(Split_Style-1)) + Split_Style^Generation_Count;

% Which nodes connect to each edge?
Node_L_Array = [1, reshape(repmat(2:(1+(Split_Style^Generation_Count - 1)/(Split_Style-1)),Split_Style,1),1,[]), (2+(Split_Style^Generation_Count - 1)/(Split_Style - 1)):(Num_Nodes - 1)];
Node_R_Array = [2:(Num_Nodes - (1 + (Split_Style^Generation_Count - 1)/(Split_Style-1))), reshape(repmat((Num_Nodes - (1 + (Split_Style^Generation_Count - 1)/(Split_Style-1)) + 1):(Num_Nodes - 1),Split_Style,1),1,[]), Num_Nodes];

% Which generation is each edge in?
Generation_Index_Array = [];
for i = 1:(Generation_Count+1)
    Generation_Index_Array(end+1:(end+Split_Style^(i-1))) = i;
end
for i = 1:(Generation_Count+1)
    Generation_Index_Array(end+1:(end+Split_Style^(Generation_Count-i+1))) = Generation_Count+i+1;
end

% Assign lengths & areas based on each edge's generation
Length_Array = generation_lengths(Generation_Index_Array);
Area_Array = generation_areas(Generation_Index_Array);

% Initialize the peaks
% any peaks must be initialized as nonzero
% any peak that hasn't completed its path currently exists
Peak_x_Array = Length_Array;
Peak_x_Array(1) = 7e-3;

Has_Peak_Array = (Peak_x_Array > 0) & (Peak_x_Array < Length_Array);
Num_Peaks_Array = Has_Peak_Array;

% Remove any peaks that don't actually exist yet
for edge = 1:Num_Edges
    source_node = Node_L_Array(edge); % the node this channel started from
    parent_edges = find(Node_R_Array == source_node); % get the parent(s) of this channel
    
    % if any parent doesn't have a peak yet or it hasn't gone the whole way, this current channel doesn't get a peak yet
    for parent = parent_edges
        if(Peak_x_Array(parent) < Length_Array(parent))
            Peak_x_Array(edge) = 0;
            break;
        end
    end
end

% Store when and where each peak was
peak_data = cell(Num_Edges, 2);
channel_ancestry_lengths = zeros(1, Num_Edges); % get the total distance the channel is from the beginning

for edge = 2:Num_Edges % skip edge 1 because it doesn't have any parents
    source_node = Node_L_Array(edge); % the node this channel started from
    parent_edge_1 = find(Node_R_Array == source_node,1); % get the first parent of this channel
    channel_ancestry_lengths(edge) = channel_ancestry_lengths(parent_edge_1) + Length_Array(parent_edge_1);
end

Plotting_Array=zeros(2, size(Peak_x_Array, 2), 2);
 
I_Array=zeros(Num_Edges,1);
V_Array=zeros(Num_Edges,1);

Channel_Array=[];
Node_Array=[];

for i=1:Num_Edges 
    Channel_Array=[Channel_Array;Channel(Length_Array(i),Area_Array(i), Node_L_Array(i), Node_R_Array(i), Has_Peak_Array(i), Num_Peaks_Array(i), Peak_x_Array(i))];
end

% compute the resistance of each channel
for i=1:Num_Edges
   Compute_R(Channel_Array(i));
end

for i=1:Num_Nodes 
    Node_Array=[Node_Array;Node()];
end

% setting the edges for each node (left and right)
for i=1:Num_Edges
    Node_Array(Channel_Array(i).Node_L).Edges_R=[Node_Array(Channel_Array(i).Node_L).Edges_R;i];
    Node_Array(Channel_Array(i).Node_L).Num_Edges_R=Node_Array(Channel_Array(i).Node_L).Num_Edges_R+1;
    Node_Array(Channel_Array(i).Node_R).Edges_L=[Node_Array(Channel_Array(i).Node_R).Edges_L;i];
    Node_Array(Channel_Array(i).Node_R).Num_Edges_L=Node_Array(Channel_Array(i).Node_R).Num_Edges_L+1;
end

[I_Array, V_Array]=Solve_Network(Channel_Array, Node_Array, Num_Edges, Num_Nodes, Network_V, Network_I);

% time stepping portion
T_final = 5E3;
dt = 10e-1;

t=0;
count=1;
time(1)=0;

count_Completed = count; % timestep that the last peak reaches the end.

while t<T_final
    % Compute resistances & solve electrical network
    for i=1:Num_Edges
        Compute_R(Channel_Array(i));
    end
    [I_Array, V_Array]=Solve_Network(Channel_Array, Node_Array, Num_Edges, Num_Nodes, Network_V, Network_I);
    
    for i=1:Num_Edges
        % Advance peaks
        if(Channel_Array(i).Has_Peak)
            if Channel_Array(i).Num_Peaks == 1
                Channel_Array(i).Peak_V = Channel_Array(i).Mu_L*I_Array(i)/(Channel_Array(i).Area*Channel_Array(i).Sigma_L);
                Channel_Array(i).Peak_x = Channel_Array(i).Peak_x+Channel_Array(i).Peak_V*dt;
            else                
                Channel_Array(i).Peak_V(1) = Channel_Array(i).Mu_L*I_Array(i)/(Channel_Array(i).Area*Channel_Array(i).Sigma_L);
                Channel_Array(i).Peak_V(2) = Channel_Array(i).Mu_L*I_Array(i)/(Channel_Array(i).Area*Channel_Array(i).Sigma_i);
                Channel_Array(i).Peak_x(1) = Channel_Array(i).Peak_x(1)+Channel_Array(i).Peak_V(1)*dt;
                Channel_Array(i).Peak_x(2) = Channel_Array(i).Peak_x(2)+Channel_Array(i).Peak_V(2)*dt;
            end 
        end
        
        % Transition peaks
        if(Channel_Array(i).Has_Peak)
            % peak merge
            if (Channel_Array(i).Num_Peaks == 2)
                if Channel_Array(i).Peak_x(2) > Channel_Array(i).Peak_x(1)
                    Channel_Array(i).Num_Peaks = Channel_Array(i).Num_Peaks - 1;
                    Channel_Array(i).Peak_x(2) = [];
                    Channel_Array(i).Peak_V(2) = [];
                end 
            end 

            % transition to next channels
            if Channel_Array(i).Peak_x(1) > Channel_Array(i).Length
                Channel_Array(i).Num_Peaks = 0;
                Channel_Array(i).Has_Peak = false;
                Channel_Array(i).Peak_x(1) = Channel_Array(i).Length;

                for j=1:Node_Array(Channel_Array(i).Node_R).Num_Edges_R
                    if ~Channel_Array(Node_Array(Channel_Array(i).Node_R).Edges_R(j)).Has_Peak
                        Channel_Array(Node_Array(Channel_Array(i).Node_R).Edges_R(j)).Has_Peak = true;
                        Channel_Array(Node_Array(Channel_Array(i).Node_R).Edges_R(j)).Peak_x(1) = 0;
                        Channel_Array(Node_Array(Channel_Array(i).Node_R).Edges_R(j)).Num_Peaks = 1;
                    else
                        Channel_Array(Node_Array(Channel_Array(i).Node_R).Edges_R(j)).Peak_x(2) = 0;
                        Channel_Array(Node_Array(Channel_Array(i).Node_R).Edges_R(j)).Num_Peaks = 2;
                    end 
                end
            end 
        end 
    end

    Peaks = [];

    for i=1:Num_Edges
       if(Channel_Array(i).Has_Peak)
           if Channel_Array(i).Num_Peaks == 1
               Peaks=[Peaks [Channel_Array(i).Peak_x(1); 0]];
           else 
               Peaks=[Peaks Channel_Array(i).Peak_x(1:2)'];
           end 

           for k=1:min(Channel_Array(i).Num_Peaks, 2)
               peak_data{i, k}(:,end+1) = [t; channel_ancestry_lengths(i) + Channel_Array(i).Peak_x(k)'];
           end 
       else
           if max(Channel_Array(i).Peak_x) > 0
               Peaks=[Peaks [Channel_Array(i).Length; Channel_Array(i).Length]];
           else
               Peaks=[Peaks [0; 0]];
           end 
       end

    end

    Plotting_Array = cat(3, Plotting_Array, Peaks);
    t=t+dt;
    count=count+1;
    time(count)=t;
    
    if (all(Peaks(1, :) >= Length_Array) && count_Completed <= 1) % If all peaks have reached the end of their channels
        count_Completed = count;
    end
end

%%

% Visualizing the whole network
channel_gap = 1;
    % gap between channels - no physical meaning, just for visualization
phi_array = Plotting_Array ./ Length_Array;
    % progress along each channel
generation_gap = channel_gap * [0, Split_Style.^((Generation_Count - 1):-1:0), -Split_Style.^(0:(Generation_Count - 1)), 0];
    % the channels in each generation will leave gaps between nodes of
    % these amounts
generation_climb_factor = -(1:Split_Style) + (Split_Style + 1)/2;
    % some channels move up, some move down, and some move more than
    % others. within a generation, the first one will move up the most, the
    % last will move down the most.
channel_climb_factor = [0, repmat(generation_climb_factor, 1, 2*((Split_Style^Generation_Count - 1)/(Split_Style - 1))), 0];
    % some channels move up, some stay 
channel_ascent = generation_gap(Generation_Index_Array) .* channel_climb_factor;
    % the distance that each channel will move up or down
channel_extend_vector = [Length_Array; channel_ascent];
    % each channel moves to the right an amount equal to length and up
    % equal to the ascent.

% now compute the point where each channel starts
channel_origin = zeros(size(channel_extend_vector));

for edge = 2:Num_Edges % skip edge 1 because it doesn't have any parents
    source_node = Node_L_Array(edge); % the node this channel started from
    parent_edge_1 = find(Node_R_Array == source_node,1); % get the first parent of this channel
    channel_origin(:,edge) = channel_origin(:,parent_edge_1) + channel_extend_vector(:,parent_edge_1);
        % each channel starts where its parent ends
end

% Now need to make time series of all edges and their points in time.
% Should be an x & y value for every edge at every point in time. So the
% resulting matrix should be 4D: <num timesteps> x <edge count> x <peak number> x 2

peak_position(:,:,1,1) = squeeze(Plotting_Array(1, :, :)) + channel_origin(1,:)'; % x(t)
peak_position(:,:,1,2) = (channel_ascent .* squeeze(phi_array(1, :, :))' + channel_origin(2,:))'; % y(t)

peak_position(:,:,2,1) = squeeze(Plotting_Array(2, :, :)) + channel_origin(1,:)'; % x(t)
peak_position(:,:,2,2) = (channel_ascent .* squeeze(phi_array(2, :, :))' + channel_origin(2,:))'; % y(t)


% Make & setup the figure
xMin = min(min(peak_position(:,:,1),[],'all'),0);
xMax = max(peak_position(:,:,1),[],'all');
yMin = min(peak_position(:,:,2),[],'all');
yMax = max(peak_position(:,:,2),[],'all');
axis([0, xMax, yMin * 1.1, yMax * 1.1]);

% leaving a little buffer above and below
set(gcf, 'Position', [100, 100, 1400, 300]);

title('Network in Time');
xlabel('Peak Position');

dStep = 20; % plot dStep timesteps at once to go faster
tMax = length(time);

for step = 1:dStep:count_Completed
    clf
    t = time(step);
    for edge = 1:Num_Edges
        plot(channel_origin(1, edge) + [0, channel_extend_vector(1, edge)], ...
            channel_origin(2, edge) + [0, channel_extend_vector(2, edge)], 'k'); hold on
        
        if ~isempty(peak_data{edge})
            stepIdx = find(peak_data{edge, 1}(1,:) == t);
            if ~isempty(stepIdx)
                plot(peak_data{edge, 1}(2,stepIdx), peak_position(edge, step, 1, 2), 'ro');
                hold on
            end

            if ~isempty(peak_data{edge, 2})
                stepIdx = find(peak_data{edge, 2}(1,:) == t);
                if ~isempty(stepIdx)
                    plot(peak_data{edge, 2}(2,stepIdx), peak_position(edge, step, 2, 2), 'bo');
                    hold on
                end
            end 
        end
    end
    title('Time: ' + string(t) + ' s')
    box off
    axis off
    drawnow
end