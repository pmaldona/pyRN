<script>
	export let text;
    export let action;
    export let style = {class: "waves-effect waves-light btn", style: "margin: 24px;"};
    export let active = true;
</script>

<!-- svelte-ignore a11y-missing-attribute -->
<!-- svelte-ignore a11y-click-events-have-key-events -->
<a 
    class={style.class + (active ? "" : " disabled")} 
    style={style.style} 
    on:click={action}
>
    {text}
</a>