<script>
    //import { ipcRenderer } from "electron";

	export let initialValues;
    export let setName;

    window.eel.expose(say_hello_js, 'say_hello_js');
	function say_hello_js(x) {
		console.log("Hello from " + x);
	}

    function openFile() {
        window.electron.open().then(result => {
            let _name = result;
            let success = eel.openFile(result)();
            success.then(result => {
                console.log(result);
                if(result == true) {
                    setName(_name);
                }
            });
            
        });
    }
</script>

<main>
	<h1>CRNS UI</h1>

    <!-- svelte-ignore a11y-missing-attribute -->
    <a class="waves-effect waves-light btn" on:click={openFile}>Open File</a>
</main>

<style>
	main {
		text-align: center;
		padding: 1em;
		max-width: 240px;
		margin: 0 auto;
	}

	h1 {
		color: #ff3e00;
		text-transform: uppercase;
		font-size: 4em;
		font-weight: 100;
	}

	@media (min-width: 640px) {
		main {
			max-width: none;
		}
	}
</style>