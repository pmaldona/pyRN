import { writable } from 'svelte/store';
import { generateGraph } from '../actions/drawNetwork';

function createSynStr() {
    const { subscribe, set, update } = writable({nodes: [], edges: [], options: {}});

    return {
		subscribe: () => subscribe(s => {
            if(s.nodes) {
                if(s.nodes.length != 0) {
                    generateGraph('hasse', s);
                }
                
            }
        }),
		set_syn_str: (syn) => set(syn),
        update_nodes: (nodes) => update(s => {
            s.nodes = nodes;
        })
	};
}

function layout(nodes, x_space = 75, y_space = 75) {
    let length_dict = {}
    nodes.forEach(node => {
        let len = node.title.split(" ").length;
        if(length_dict[len]){
            length_dict[len]["length"] = length_dict[len]["length"] + 1;
        }
        else{
            length_dict[len] = {};
            length_dict[len]["c"] = 0;
            length_dict[len]["length"] = 1;
        }
    });

    console.log(length_dict);

    nodes.forEach(node => {
        console.log(node.title + " => " + node.title.split(" ").length);
        let len = node.title.split(" ").length;
        let count = length_dict[len]["length"];
        let _y = -y_space * (len-1);
        let _x = 0;
        if (count > 1){
            _x = -x_space*((count-1)/2) + length_dict[len]["c"]*x_space
            length_dict[len]["c"] = length_dict[len]["c"]+1;
        }
        
        node["y"] = _y;
        node["x"] = _x;
        node["fixed"] = {x: false, y: true}
        console.log(length_dict);
    });

    return nodes;
}

export const synStr = createSynStr();
export { layout };