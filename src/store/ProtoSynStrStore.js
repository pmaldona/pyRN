import { writable } from 'svelte/store';
import { generateGraph } from './../actions/drawNetwork';

function createProtoSyn() {
    const { subscribe, set, update } = writable({nodes: [], edges: [], options: {}});

    return {
		subscribe: () => subscribe(p => {
            if(p.nodes) {
                if(p.nodes.length != 0) {
                    generateGraph('network', p);
                }
                
            }
            console.log(p);
        }),
		set_proto_syn: (proto) => set(proto),
        update_nodes: (nodes) => update(p => {
            console.log(nodes);
            console.log(p);
            p.nodes = nodes;
        })
	};
}

export const protoSyn = createProtoSyn();