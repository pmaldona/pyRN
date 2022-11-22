<script>
	import Base from '../misc/Base.svelte';
  	import Network from '../Network/Network.svelte';
	import Hasse from '../Hasse/Hasse.svelte';
	import Simulation from '../Simulation/Simulation.svelte';

  	const pages = [Base, Network, Hasse, Simulation];

	let page = 0;

	const is_loaded = async () => {
        let promise = eel.has_loaded_file()();
        let bool = false;
		await promise.then(result => {
            bool = result;
        });
		console.log(bool);
		return bool;
    }

	const openPage = (_p) => {
		page = _p;
		console.log(page);
	}
</script>

<nav>
	<div class="nav-wrapper">
		<ul id="nav-mobile" class="left">
			<!-- svelte-ignore a11y-missing-attribute -->
			<li on:click={() => {openPage(0)}}><a>Main</a></li>
			<!-- svelte-ignore a11y-missing-attribute -->
			<li on:click={() => {openPage(1)}}><a>Reaction Network</a></li>
			<!-- svelte-ignore a11y-missing-attribute -->
			<li on:click={() => {openPage(2)}}><a>Hasse Diagramm</a></li>
			<!-- svelte-ignore a11y-missing-attribute -->
			<li on:click={() => {openPage(3)}}><a>Simulation</a></li>
		</ul>
	</div>
</nav>

<svelte:component
  	this={pages[page]}
	openPage={openPage}
	is_loaded={is_loaded}
/>

<!-- <footer>
	<p>{ $filename || "No file opened."}</p>
</footer> -->
<!-- 
<style>
	footer {
		position: absolute;
		bottom: 0px;
	}
</style> -->