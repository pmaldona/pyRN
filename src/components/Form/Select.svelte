<script>
    export let bind;
    export let callback;
    export let elements;
    export let style = "display:block;";
</script>

<select bind:value={bind} on:change="{() => callback}" style={style}>
    {#each elements as element}
        <option value={element}>
            {element}
        </option>
    {/each}
</select>