using Graphs

# find the connected graph containing the slack bus
function find_connected_graph(nodes, edges)
	n = length(nodes)
	A = zeros(Int64, n,n)
    	for edge in edges
			if edge.line_status
				A[edge.source.id, edge.target.id] = 1
				A[edge.target.id, edge.source.id] = 1
			end
    	end
	queue = filter(n -> n.bus_type == 3, nodes)[1].id
	connected = [];
	while(length(queue) != 0)
        connected  = [connected; queue[1]];
		x = findin(A[:,queue[1]],1)
        queue = [queue; setdiff(x, connected)];
        queue = unique(queue[2:end]);
    end
	#is_connected = falses(n,1)
	#is_connected[connected] = true
	return sort(connected)
end


