import {writable} from "svelte/store";

function activate(initState){
    let {subscribe, update} = writable(initState);
    const activate = () => update(s => true)
    return {subscribe, activate}
}

export const createActivation = (state) => activate(state)