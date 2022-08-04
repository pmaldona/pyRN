import { writable } from 'svelte/store';
import { generateGraph } from './../actions/drawNetwork';

function createNetwork() {
    const { subscribe, set, update } = writable({nodes: [], edges: [], options: {}});

    return {
		subscribe: () => subscribe(n => {
            if(n.nodes) {
                if(n.nodes.length != 0) {
                    generateGraph('network', n);
                }
                
            }
            console.log(n);
        }),
		set_network: (ntwrk) => set(ntwrk),
        update_nodes: (nodes) => update(n => {
            console.log(nodes);
            console.log(n);
            n.nodes = nodes;
        })
	};
}

export const network = createNetwork();