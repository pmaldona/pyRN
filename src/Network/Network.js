import { setVisObject } from '../misc/drawNetwork';

export function genRNStr(type) {
    let promise = eel.gen_network()();
    let result = promise.then(result => {
        if (result == null) {
            return;
        }
        for(let i = 0; i < result.edges.length; i++) {
            let edge = result.edges[i];
            if(edge.title){
                edge.label = edge.title;
            }
            edge.color = {color: '#848484', highlight: '#848484', hover: '#848484', inherit: false, opacity: 1.0};
            edge.smooth = {enabled: false};
            for(let j = 0; j < result.edges.length; j++) {
                if(j != i) {
                    if(edge.from == result.edges[j].to && edge.to == result.edges[j].from) {
                        edge.smooth = {enabled: true, type: "discrete", roundness: 0.7};
                    }
                }
            }
            edge.physics = false;
        }
        // if(type == 'org_network') {
        //     console.log(result);
        //     return result;
        // }
        return setVisObject(type, result.nodes, result.edges, result.options, true);
    }); 
    return result;
}