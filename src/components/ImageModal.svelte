<script context="module" lang="ts">
	let onTop   //keeping track of which open modal is on top
	const modals={}  //all modals get registered here for easy future access
	
	// 	returns an object for the modal specified by `id`, which contains the API functions (`open` and `close` )
	export function getModal(id=''){
		return modals[id]
	}
</script>

<script lang="ts">
	import {onDestroy} from 'svelte'
	
	let topDiv
	let visible=false
	let prevOnTop
	let closeCallback

    let image_source;

	export let id=''
	export let heading="";
	export let contentStyle="";

    function keyPress(ev){
        //only respond if the current modal is the top one
        if(ev.key=="Escape" && onTop==topDiv) close(0) //ESC
    }

    /**  API **/
    function open(callback){
        closeCallback=callback
        if(visible) return
        prevOnTop=onTop
        onTop=topDiv
        window.addEventListener("keydown",keyPress)
        
        //this prevents scrolling of the main window on larger screens
        document.body.style.overflow="hidden" 

        visible=true
        //Move the modal in the DOM to be the last child of <BODY> so that it can be on top of everything
        document.body.appendChild(topDiv)
    }
        
    function close(retVal){
        if(!visible) return
        window.removeEventListener("keydown",keyPress)
        onTop=prevOnTop
        if(onTop==null) document.body.style.overflow=""
        visible=false
        if(closeCallback) closeCallback(retVal)
    }
        
    //expose the API
    modals[id]={open,close}
        
    onDestroy(()=>{
        delete modals[id]
        window.removeEventListener("keydown",keyPress)
    })

    export async function init(img) {
        image_source = img;
    }
	
</script>

<main>
    {#if image_source != undefined}
        <!-- svelte-ignore a11y-missing-attribute -->
        <img src='data:image/png;base64,{image_source}' style="max-width: 95%; height: auto; object-fit: contain;">
    {/if}
</main>

<style>
	#topModal {
		visibility: hidden;
		z-index: 9999;
		position: fixed;
		top: 0;
		left: 0;
		right: 0;
		bottom: 0;
		background: #4448;
		display: flex;
		align-items: center;
		justify-content: center;
	}
	#modal {
		position: relative;
		border-radius: 4px;
		background: white;
    border: 2px solid rgb(40, 40, 40)0;
		/* filter: drop-shadow(5px 5px 5px #555); */
		padding: 1em;
	}

	.visible {
		visibility: visible !important;
	}

	#close {
		position: absolute;
		top:-12px;
		right:-12px;
		width:24px;
		height:24px;
		cursor: pointer;
		fill:#F44;
		transition: transform 0.3s;
	}	

	#close:hover {
		transform: scale(2);
	}

	#close line {
		stroke:#FFF;
		stroke-width:2;
	}
	#modal-content {
		max-width: calc(100vw - 20px);
		max-height: calc(100vh - 20px);
		overflow: auto;
	}
</style>