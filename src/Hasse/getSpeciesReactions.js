import { GraphType } from '../misc/drawNetwork';
import { genRNStr } from '../Network/Network';

export default async function getSpeciesReactions(organization) {
    console.log(organization.title);
    let species_ids = organization.title.substring(1,organization.title.length-1).replaceAll("'", "").replaceAll("\n","").replaceAll("\r","").split(" ");
    let reaction_ids = ["id1"];
    let reaction_str = await eel.get_reactions()();
    let reactions = [];

    let network = await genRNStr(GraphType.OrgNetwork);
    let nodes = Object.values(network.body.nodes);

    nodes.forEach(node => {
        if(!species_ids.includes(node.id)){
            let included = true
            if(node.options.shape=="square") {
                node.edges.forEach(edge => {
                    let id = edge.toId;
                    if(node.id == edge.toId) {
                        id = edge.fromId;
                    }

                    if(species_ids.includes(id) == false) {
                        included = false;
                    }
                });
            }else {
                included = false;
            }
            if(included == false) {
                
            }
            else {
                if(node.options.shape=="square") {
                    if(!reaction_ids.includes(node.id)) {
                        let reac = "";
                        reaction_str.forEach(str => {
                            if (str.substring(0,node.id.length)== node.id) {
                                reac = str.substring(node.id.length+2);
                            }
                        });
                        // included_reactions.push(node);
                        reactions.push({id: node.id, str: reac});
                    }
                }
            }
        }
    });

    return {species: species_ids, reactions: reactions};
}